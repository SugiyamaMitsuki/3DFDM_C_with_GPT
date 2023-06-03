#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct ReceiverInfo {
    int id;
    double lat;
    double lon;
    double depth;// m 
    std::string name;
};

std::vector<ReceiverInfo> readReceiversCSV(const std::string& filename) {
    std::vector<ReceiverInfo> receivers;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '#') {
            continue;
        }
        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        std::istringstream ss(line);
        ReceiverInfo receiver;
        char separator;

        ss >> receiver.id >> separator >> receiver.lat >> separator >> receiver.lon >> separator >> receiver.depth >> separator >> receiver.name;
        receiver.depth = receiver.depth * 1000.0;

        receivers.push_back(receiver);
    }

    return receivers;
}

