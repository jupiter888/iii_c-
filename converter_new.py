import sys
import os
import dask.dataframe as dd
import pandas as pd
import csv
import pyarrow as pa
import pyarrow.parquet as pq

def detect_delimiter(filepath, sample_size=1000):
    with open(filepath, 'r', encoding='utf-8') as f:
        sample = [next(f) for _ in range(sample_size)]
    sniffer = csv.Sniffer()
    return sniffer.sniff("\n".join(sample)).delimiter

def convert_csv_to_parquet(csv_path, output_dir):
    print(f"🔍 validating csv before conversion: {csv_path}")

    # detect delimiter
    delimiter = detect_delimiter(csv_path)
    print(f"🧐 detected delimiter: {repr(delimiter)}")
    
    # ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # specify columns to keep and enforce dtypes
    keep_columns = [
        'gbifID', 'species', 'decimalLongitude', 'decimalLatitude',
        'countryCode', 'elevation', 'datasetKey', 'eventDate'
    ]
    dtype_fix = {
        'gbifID': 'object',
        'species': 'object',
        'decimalLongitude': 'object',  # will convert later
        'decimalLatitude': 'object',   # will convert later
        'countryCode': 'object',
        'elevation': 'object',         # will convert later
        'datasetKey': 'object',
        'eventDate': 'object'
    }
    
    try:
        ddf = dd.read_csv(
            csv_path,
            delimiter=delimiter,
            usecols=keep_columns,
            dtype=dtype_fix,
            assume_missing=True
        )
    except Exception as e:
        print(f"❌ error reading csv: {e}")
        sys.exit(1)
    
    # convert specific columns to numeric
    for col in ['decimalLongitude', 'decimalLatitude', 'elevation']:
        ddf[col] = ddf[col].apply(pd.to_numeric, errors='coerce', meta=(col, 'float64'))
    
    # check partition consistency
    print("\n📊 partition report:")
    partition_info = ddf.map_partitions(lambda df: pd.DataFrame({'row_count': [len(df)]})).compute()
    print("Row counts per partition:")
    print(partition_info)
    print("Dtypes:")
    print(ddf.dtypes)
    
    # missing value report
    missing_report = ddf.isnull().sum().compute()
    print("\n📉 missing values per column before conversion:")
    print(missing_report[missing_report > 0])
    
    # fill missing values
    ddf = ddf.fillna({
        'species': 'Unknown',
        'countryCode': 'Unknown',
        'decimalLongitude': -9999,
        'decimalLatitude': -9999,
        'elevation': -9999,
        'eventDate': 'Unknown'
    })
    
    try:
        ddf.to_parquet(output_dir, engine='pyarrow', write_index=False)
        print(f"✅ conversion complete. Parquet files saved in: {output_dir}")
    except Exception as e:
        print(f"❌ error writing Parquet files: {e}")
        sys.exit(1)
    
    # generate a report for each parquet file
    print("\n📑 parquet file report:")
    for filename in os.listdir(output_dir):
        if filename.endswith(".parquet"):
            file_path = os.path.join(output_dir, filename)
            try:
                table = pq.read_table(file_path)
                num_rows = table.num_rows
                schema = table.schema
                print(f"File: {filename} - Rows: {num_rows}")
                print("Schema:")
                print(schema)
            except Exception as e:
                print(f"❌ error reading parquet file {filename}: {e}")
    
    # validate each parquet file by attempting to open it
    print("\n🔎 validating each parquet file:")
    for filename in os.listdir(output_dir):
        if filename.endswith(".parquet"):
            file_path = os.path.join(output_dir, filename)
            try:
                _ = pq.read_table(file_path)
                print(f"File {filename} is valid.")
            except Exception as e:
                print(f"❌ file {filename} is invalid: {e}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python converter_newer.py <input_csv_file> <output_parquet_dir>")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_dir = sys.argv[2]

    convert_csv_to_parquet(csv_file, output_dir)
