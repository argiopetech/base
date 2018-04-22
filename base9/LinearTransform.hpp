#ifndef LINEAR_TRANSFORM_HPP
#define LINEAR_TRANSFORM_HPP

#include <cassert>
#include <utility>

#include <stdexcept>

#include <iostream>
using std::cout;
using std::endl;

enum class Direction { None, Low, High };
enum class TransformMethod { Interp, Extrap, ExtrapLow, ExtrapHigh };

// This needed to be rather flipped inside-out to work with C++'s view of objects.
// data Result = OnGrid Double
//             | Pinned Direction Double
//             | Extrapolated Direction Double
template <TransformMethod M>
struct Result
{
    constexpr Result() = delete;
};

template <>
struct Result<TransformMethod::Interp>
{
    Result(double newVal, bool, bool = false)
        : val(newVal)
    {
        throw std::logic_error("This constructor should never be called.");
    }

    constexpr Result(double newVal)
        : val(newVal)
    {}

    constexpr Result(double newVal, Direction newDir)
        : isPinned(true), dir(newDir), val(newVal)
    {}

    const bool isPinned = false;
    const Direction dir = Direction::None;

    const double val;
};

template <>
struct Result<TransformMethod::Extrap>
{
    Result(double newVal, bool, bool = false)
        : val(newVal)
    {
        throw std::logic_error("This constructor should never be called.");
    }

    constexpr Result(double newVal)
        : val(newVal)
    {}

    constexpr Result(double newVal, Direction newDir)
        : isExtrapolated(true), dir(newDir), val(newVal)
    {}

    const bool isExtrapolated = false;
    const Direction dir = Direction::None;

    const double val;
};

template <>
struct Result<TransformMethod::ExtrapLow>
{
    Result(double newVal, Direction)
        : val(newVal)
    {
        throw std::logic_error("This constructor should never be called");
    }

    constexpr Result(double newVal, bool pinned = false, bool extrapped = false)
        : isPinned(pinned), isExtrapolated(extrapped), val(newVal)
    {}

    const bool isPinned = false;
    const bool isExtrapolated = false;

    const double val;
};

template <>
struct Result<TransformMethod::ExtrapHigh>
{
    Result(double newVal, Direction)
        : val(newVal)
    {
        throw std::logic_error("This constructor should never be called.");
    }

    constexpr Result(double newVal, bool pinned = false, bool extrapped = false)
        : isPinned(pinned), isExtrapolated(extrapped), val(newVal)
    {}

    const bool isPinned = false;
    const bool isExtrapolated = false;

    const double val;
};

// Template metaprogramming. The TransformMethod, being known at runtime,
// with a default value of "Extrap"
template <TransformMethod M = TransformMethod::Extrap>
auto linearTransform(double xL, double xH, double yL, double yH, double xActual) -> const Result<M>
{
    // Defer calculation till we know we need it
    // Shares function parameters without being visible externally.
    auto linInterpExtrap = [=]()
    {
        double position = (xActual - xL) / (xH - xL);
        return yL + position * (yH - yL);
    };

    // Could have done this in a switch statement, but if's seem cleaner for the long term case
    if (M == TransformMethod::Interp)
    {
        // Interp | xActual > xH -> Pinned High yH
        //        | xActual < xL -> Pinned Low  yL
        //        | otherwise    -> OnGrid linInterpExtrap
        if (xActual > xH)
            return {yH, Direction::High};
        else if (xActual < xL)
            return {yL, Direction::Low};
        else
            return linInterpExtrap();
    }
    else if (M == TransformMethod::Extrap)
    {
        // Extrap | xActual > xH -> Extrapolated High linInterpExtrap
        //        | xActual < xL -> Extrapolated Low  linInterpExtrap
        //        | otherwise    -> OnGrid linInterpExtrap
        if (xActual > xH)
            return {linInterpExtrap(), Direction::High};
        else if (xActual < xL)
            return {linInterpExtrap(), Direction::Low};
        else
            return linInterpExtrap();
    }
    else if (M == TransformMethod::ExtrapLow)
    {
        // (ExtrapOnly Low)  | xActual > xH -> Pinned High yH
        //                   | xActual < xL -> Extrapolated Low linInterpExtrap
        //                   | otherwise    -> OnGrid linInterpExtrap
        if (xActual > xH)
            return {yH, true};
        else if (xActual < xL)
            return {linInterpExtrap(), false, true};
        else
            return linInterpExtrap();
    }
    else if (M == TransformMethod::ExtrapHigh)
    {
        // (ExtrapOnly High) | xActual > xH -> Extrapolated High linInterpExtrap
        //                   | xActual < xL -> Pinned Low yL
        //                   | otherwise    -> OnGrid linInterpExtrap
        if (xActual > xH)
            return {linInterpExtrap(), false, true};
        else if (xActual < xL)
            return {yL, true};
        else
            return linInterpExtrap();
    }
    else
    {
        throw std::logic_error("Invalid TransformMethod in linearTransform.");
    }
}

template <TransformMethod M>
auto linearTransform(const std::pair<double, double> &x, const std::pair<double, double> &y, double p) -> Result<M>
{
    return linearTransform<M>(x.first, x.second, y.first, y.second, p);
}

template <TransformMethod M, class ForwardIterator1, class ForwardIterator2>
auto linearTransform(const ForwardIterator1 &it1, const ForwardIterator2 &it2, double p) -> Result<M>
{
    return linearTransform<M>(it1[0], it1[1], it2[0], it1[1], p);
}

#endif
