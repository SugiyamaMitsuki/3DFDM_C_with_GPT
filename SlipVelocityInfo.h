#ifndef COLUMNAR_SLIP_READER_H
#define COLUMNAR_SLIP_READER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

class SlipVelocityReader
{
public:
    struct SlipVelocity
    {
        std::vector<double> time;
        std::vector<double> velocity;
    };

    explicit SlipVelocityReader(const std::string &filename)
    {
        readCSV(filename);
    }

    const std::vector<SlipVelocity> &getData() const
    {
        return data;
    }

private:
    std::vector<SlipVelocity> data;

    void readCSV(const std::string &filename)
    {
        std::ifstream file(filename);

        if (!file.is_open())
        {
            throw std::runtime_error("Could not open the file: " + filename);
        }

        std::string line;
        // Skip header lines
        for (int i = 0; i < 3; ++i)
        {
            std::getline(file, line);
        }

        std::vector<std::vector<double>> columns;

        while (std::getline(file, line))
        {
            if (line.empty())
            {
                continue;
            }

            std::istringstream ss(line);
            double value;
            char separator;

            size_t columnIndex = 0;
            while (!ss.eof())
            {
                ss >> value >> separator;
                if (columns.size() <= columnIndex)
                {
                    columns.emplace_back();
                }
                columns[columnIndex].push_back(value);
                ++columnIndex;

                // //debug
                // std::cout << value << " ";
                // std::cout << separator << " ";
                // std::cout << columnIndex << " ";
                // std::cout << std::endl;
            }
        }

        SlipVelocity slipVelocity;
        slipVelocity.time = columns[0];

        for (size_t i = 1; i < columns.size(); ++i)
        {
            slipVelocity.velocity.clear();

            for (size_t j = 0; j < columns[i].size(); ++j)
            {
                slipVelocity.velocity.push_back(columns[i][j]);
            }
            data.push_back(slipVelocity);
        }
    }
};

#endif // COLUMNAR_SLIP_READER_H
