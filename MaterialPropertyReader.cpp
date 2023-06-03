#include "MaterialPropertyReader.h"
#include <iostream>
#include <fstream>
#include <sstream>

MaterialPropertyReader::MaterialPropertyReader(const std::string &filename) : file(filename)
{
    loadData();
}

void MaterialPropertyReader::loadData()
{
    std::ifstream inputStream(file);
    std::string line;

    if (inputStream.is_open())
    {
        // Skip header lines
        int headerLines = 0; // 7;
        for (int i = 0; i < headerLines; ++i)
        {
            std::getline(inputStream, line);
        }

        // Read data
        while (std::getline(inputStream, line))
        {
            std::stringstream lineStream(line);
            std::string cell;
            MaterialProperty row;

            // 1文字目が#の行はコメント行として無視する
            if (line[0] == '#')
            {
                continue;
            }
            // Skip empty lines
            if (line.empty())
            {
                continue;
            }
            // Read data
            std::getline(lineStream, cell, ',');
            row.STN = std::stoi(cell);
            std::getline(lineStream, cell, ',');
            row.SVP = std::stod(cell);
            std::getline(lineStream, cell, ',');
            row.SVS = std::stod(cell);
            std::getline(lineStream, cell, ',');
            row.SRO = std::stod(cell);
            std::getline(lineStream, cell, ',');
            row.SQP = std::stod(cell);
            std::getline(lineStream, cell, ',');
            row.SQS = std::stod(cell);

            // debug
            // std::cout << row.STN << ", " << row.SVP << ", " << row.SVS << ", " << row.SRO << ", " << row.SQP << ", " << row.SQS << std::endl;
            //

            data.push_back(row);
        }

        inputStream.close();
    }
    else
    {
        std::cerr << "Unable to open file: " << file << std::endl;
    }
}

MaterialProperty MaterialPropertyReader::operator[](int stn)
{
    MaterialProperty result;
    for (const auto &row : data)
    {
        if (row.STN == stn)
        {
            result = row;
            break;
        }
    }

    return result;
}
