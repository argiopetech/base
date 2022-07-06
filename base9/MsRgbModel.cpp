#include <iostream>

#include "MsRgbModel.hpp"

using std::cout;
using std::endl;
using std::string;

unsigned int MsBoundsError::counts[4] = {0, 0, 0, 0};

void MsBoundsError::countSummary()
{
    string binder = " errors with ";
    string labels[4] = {"Fe/H too high", "Fe/H too low", "age too high", "age too low"};

    for (size_t i = 0; i < 4; ++i)
    {
        if (counts[i] > 0)
        {
            cout << counts[i] << binder << labels[i] << endl;
        }
    }
}
