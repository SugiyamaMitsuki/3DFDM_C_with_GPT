#include "Config.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

Config::Config(const std::string &filename)
{
    loadCSV(filename);
}

Config::~Config()
{
}

template <>
int Config::get<int>(const std::string &key) const
{
    auto it = parameters.find(key);
    if (it == parameters.end())
        throw std::runtime_error("Key not found: " + key);
    return std::stoi(it->second);
}

template <>
double Config::get<double>(const std::string &key) const
{
    auto it = parameters.find(key);
    if (it == parameters.end())
        throw std::runtime_error("Key not found: " + key);
    return std::stod(it->second);
}

template <>
std::string Config::get<std::string>(const std::string &key) const
{
    auto it = parameters.find(key);
    if (it == parameters.end())
        throw std::runtime_error("Key not found: " + key);
    return it->second;
}

void Config::loadCSV(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file: " + filename);

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, ',') && std::getline(iss, value, ','))
            parameters[key] = value;
    }

    file.close();
}