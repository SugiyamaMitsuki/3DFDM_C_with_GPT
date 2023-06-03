
#pragma once
#include <string>
#include <map>

class Config
{
public:
    Config(const std::string &filename);
    ~Config();

    template <typename T>
    T get(const std::string &key) const;

private:
    std::map<std::string, std::string> parameters;
    void loadCSV(const std::string &filename);
};
