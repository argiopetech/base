#include <stdexcept>
#include <string>

class WDBoundsError : public std::range_error
{
//    using std::range_error::range_error;

  public:
    explicit WDBoundsError (const std::string& what_arg)
        : std::range_error(what_arg) {}

    explicit WDBoundsError (const char* what_arg)
        : std::range_error(what_arg) {}
};
