#include "MpiMcmcApplication.hpp"

int main (int argc, char *argv[])
{
    MpiMcmcApplication master;

    return master.run(argc, argv);
}
