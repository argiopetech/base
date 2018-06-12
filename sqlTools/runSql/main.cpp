#include <iostream>

#include "SqlQueryRunner.hpp"

using std::cerr;
using std::endl;


int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " TARGET.base9 \"SQL CODE\"" << endl;

        return 1;
    }

    auto sql = argv[2];

    SqlQueryRunner query(argv[1], sql);

    query.run();
    query.prettyPrint();

    return 0;
}
