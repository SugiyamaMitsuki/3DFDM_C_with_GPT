#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

struct SourceInfo {
    int id;
    double M0;
    double strike;
    double dip;
    double rake;
    double magnitude;
    double lat;
    double lon;
    double depth;// m 
};

std::vector<SourceInfo> readSourcesCSV(const std::string& filename) {
    std::vector<SourceInfo> sources;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comment lines
        if (line[0] == '#') {
            continue;
        }
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        std::istringstream ss(line);
        SourceInfo src;
        char separator;
        ss >> src.id >> separator >> src.M0 >> separator >> src.strike >> separator
           >> src.dip >> separator >> src.rake >> separator >> src.magnitude >> separator
           >> src.lat >> separator >> src.lon >> separator >> src.depth;
        src.depth = src.depth * 1000.0;
        sources.push_back(src);
    }

    return sources;
}
