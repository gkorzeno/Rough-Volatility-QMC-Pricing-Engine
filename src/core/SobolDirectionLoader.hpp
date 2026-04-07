#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

struct SobolDirectionData {
    int s;
    int a;
    std::vector<int> m;
};

class SobolDirectionLoader {
public:

    static std::vector<SobolDirectionData> load(const std::string& filename)
    {
        std::ifstream file(filename);

        if(!file)
            throw std::runtime_error("Could not open Sobol direction file");

        std::vector<SobolDirectionData> table;
        std::string line;

        while(std::getline(file, line))
        {
            if(line.empty())
                continue;

            std::stringstream ss(line);

            int d, s, a;

            // skip header or malformed rows
            if(!(ss >> d >> s >> a))
                continue;

            SobolDirectionData row;
            row.s = s;
            row.a = a;

            for(int i = 0; i < s; i++)
            {
                int mi;
                ss >> mi;
                row.m.push_back(mi);
            }

            table.push_back(row);
        }

        return table;
    }
};