#include <stdexcept>

#include "MpiMcmcApplication.hpp"

using std::cerr;
using std::endl;
using std::exception;

int main (int argc, char *argv[])
{
    Settings settings;

    // Setup settings
    settings.fromCLI (argc, argv);
    if (!settings.files.config.empty())
    {
        settings.fromYaml (settings.files.config);
    }
    else
    {
        settings.fromYaml ("base9.yaml");
    }

    settings.fromCLI (argc, argv);

    MpiMcmcApplication master(settings);

    try
    {
        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\n\n" << e.what() << endl;
        return -1;
    }
}
