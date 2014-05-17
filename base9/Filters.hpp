#ifndef FILTERSET_HPP
#define FILTERSET_HPP

#include <string>
#include <unordered_map>
#include <vector>

struct Filters
{
    Filters() = delete;
    ~Filters() = default;

    static std::vector<double> calcAbsCoeffs(const std::vector<std::string>&);
    const static std::unordered_map<std::string, double> absCoeffs;
};

#endif
