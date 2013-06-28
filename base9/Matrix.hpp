#include <array>

template <class T, unsigned I, unsigned J>
using Matrix = std::array<std::array<T, J>, I>;
