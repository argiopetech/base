#include <stdexcept>

class WDBoundsError : public std::range_error
{
    using std::range_error::range_error;
};
