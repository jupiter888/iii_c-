#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <utility>

// function to load country/province centroids from a csv file
std::vector<std::pair<double,double>> loadCentroidsCSV(const std::string& filename) {
    std::vector<std::pair<double,double>> centroids;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "error: cannot open centroids csv file " << filename << std::endl;
        return centroids;
    }

    // read header
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "error: empty centroids file or missing header." << std::endl;
        return centroids;
    }

    // parse data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        double lonVal = std::numeric_limits<double>::quiet_NaN();
        double latVal = std::numeric_limits<double>::quiet_NaN();

        // example usage: read columns by name if needed,  or assume fixed column positions
        // for demonstration, let's assume columns in order: 
        // type, iso3, area, centroid.lon, centroid.lat, ...
        // we only care about centroid.lon (4th col) and centroid.lat (5th col)

        // 1. type
        if (!std::getline(ss, token, ',')) continue;
        // 2. iso3
        if (!std::getline(ss, token, ',')) continue;
        // 3. area
        if (!std::getline(ss, token, ',')) continue;
        // 4. centroid.lon
        if (!std::getline(ss, token, ',')) continue;
        try {
            lonVal = std::stod(token);
        } catch(...) {}
        // 5. centroid.lat
        if (!std::getline(ss, token, ',')) continue;
        try {
            latVal = std::stod(token);
        } catch(...) {}

        // store in vector
        if (!std::isnan(lonVal) && !std::isnan(latVal)) {
            centroids.push_back(std::make_pair(lonVal, latVal));
        }
    }

    std::cout << "loaded " << centroids.size() << " centroids from " << filename << std::endl;
    return centroids;
}
