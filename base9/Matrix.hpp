#include <array>
#include <vector>

template <class T, unsigned I, unsigned J>
using Matrix = std::array<std::array<T, J>, I>;

template <class T>
using Vatrix = std::vector<std::vector<T>>;

template <class T, unsigned I>
using MVatrix = std::array<std::vector<T>, I>;

template <class T, unsigned I>
using VMatrix = std::vector<std::array<T, I>>;
