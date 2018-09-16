#include <iostream>

#include "SqlQueryRunner.hpp"

using std::cerr;
using std::endl;


int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cerr << "Usage: " << argv[0] << " TARGET.base9" << endl;
        return 1;
    }

    auto sql =
        "select \"Multi-Pop\" as Program, runId as 'Run Identifier', time as Time, count(iterId) as Iterations"
        "  from run"
        "  join multi_pop"
        "    on run.id = runId"
        "  group by runId "

        "union "

        "select \"Single Pop\" as program, runId, time, count(iterId)"
        "  from run"
        "  join single_pop"
        "    on run.id = runId"
        "  group by runId "

        "union "

        "select \"Sample Mass\" as program, runId, time, count(distinct iterId)"
        "  from run"
        "  join sample_mass"
        "    on run.id = runId"
        "  group by runId "
        "order by runId;";

    SqlQueryRunner query(argv[1], sql);

    query.run();
    query.prettyPrint();

    return 0;
}
