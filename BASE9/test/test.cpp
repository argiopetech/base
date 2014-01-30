#include <array>
#include <iostream>
#include <string>

#include "constants.hpp"

using std::array;
using std::string;

void initBase(string, int, int, int, int);
void setClusterParameters(double, double, double, double, double, double);
void changeModels(int, int, int, int);
void setIFMRParameters(double, double, double);
array<double, FILTS> evolve (double);

int main()
{
    initBase("models/", 0, 2, 1, 1);
    setClusterParameters(8.796, 0.07, 0.0, 0.09, 0.29, 0.5);

    for (int k = 0; k < 801; ++k)
    {
        auto a = evolve(k * 0.01);
        std::cout << k * 0.01 << ":   "; 
        for (int i = 0; i < 7; ++i)
            std::cout << a.at(i) << ", ";
        std::cout << a.at(7) << std::endl;
    }

    return 0;
}
