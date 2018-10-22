template <typename T>
SqlBackingStore<T>::SqlBackingStore(const SqlBackingStore &other)
    : db(other.db), run(other.run)
{ ; }


template <typename T>
SqlBackingStore<T>::SqlBackingStore(const RunData &bootstrap)
    : db(bootstrap.db), run(bootstrap.run)
{ ; }


template <typename T>
SqlBackingStore<T>::SqlBackingStore(const std::string baseName)
{
    string dbName = baseName + ".base9";

    openDb(dbName);

    ensureTables();

    generateRun();

    buildInsertStatement();
}


template <typename T>
SqlBackingStore<T>::~SqlBackingStore()
{
    sqlite3_finalize(insertIteration);
    insertIteration = nullptr;
}


template <typename T>
void SqlBackingStore<T>::openDb(const std::string dbName)
{
    sqlite3 *tempDb;

    auto ret = sqlite3_open(dbName.c_str(), &tempDb);

    dbErrorIf(ret, "Opening database \"" + dbName + "\"");

    db.reset(tempDb, sqlite3_close);
}

template <typename T>
void SqlBackingStore<T>::ensureTables()
{
    enableForeignKeys();

    execOnlyInTransaction(
        "creating tables",
        "create table if not exists run"
        "( id   integer not null"
        ", time datetime default CURRENT_TIMESTAMP"
        ", primary key (id));"

        "create table if not exists iteration"
        "( runId integer not null"
        ", id    integer not null"
        ", primary key (runId, id)"
        ", foreign key (runId) references run);"

        "create table if not exists single_pop"
        "( runId          integer not null"
        ", iterId         integer not null"
        ", logAge         double"
        ", y              double"
        ", feh            double  not null"
        ", modulus        double  not null"
        ", absorption     double  not null"
        ", carbonFraction double"
        ", ifmrIntercept  double"
        ", ifmrSlope      double"
        ", ifmrQuadCoef   double"
        ", logPost        double  not null"
        ", stage          integer not null"
        ", primary key (runId, iterId)"
        ", foreign key (runId, iterId) references iteration);"

        "create table if not exists multi_pop"
        "( runId          integer not null"
        ", iterId         integer not null"
        ", logAge         double"
        ", feh            double  not null"
        ", modulus        double  not null"
        ", absorption     double  not null"
        ", carbonFraction double"
        ", ya             double  not null"
        ", yb             double  not null"
        ", lambda         double  not null"
        ", logPost        double  not null"
        ", stage          integer not null"
        ", primary key (runId, iterId)"
        ", foreign key (runId, iterId) references iteration);"

        "create table if not exists field_star_likelihood"
        "( runId integer not null"
        ", likelihood double not null"
        ", primary key (runId)"
        ", foreign key (runId) references run(id));"

        "create table if not exists star"
        "( runId integer not null"
        ", starId text not null"
        ", primaryMass double not null"
        ", secondaryMassRatio double not null"
        ", stage integer not null"
        ", clusterMembershipPrior double not null"
        ", useDuringBurnin boolean not null"
        ", primary key (runId, starId)"
        ", foreign key (runId) references run);"

        "create table if not exists photometry"
        "( runId integer not null"
        ", starId text not null"
        ", filter text not null"
        ", magnitude double not null"
        ", stdDev double not null"
        ", foreign key (runId, starId) references star);"

        "create table if not exists star_posterior"
        "( runId integer not null"
        ", iterId integer not null"
        ", starId text not null"
        ", clust_a_posterior double not null"
        ", clust_b_posterior double not null"
        ", primary key (runId, iterId, starId)"
        ", foreign key (runId, iterId) references iteration (runId, id)"
        ", foreign key (runId, starId) references star (runId, starId));"

        "create table if not exists sample_mass"
        "( runId integer not null"
        ", referencedRunId integer not null"
        ", iterId integer not null"
        ", starId text not null"
        ", primaryMass double not null"
        ", massRatio double not null"
        ", clusterMembership double not null"
        ", primary key (runId, referencedRunId, iterId, starId)"
        ", foreign key (referencedRunId, iterId) references iteration (runId, id)"
        ", foreign key (referencedRunId, starId) references star (runId, starId)"
        ", foreign key (referencedRunId, iterId) references single_pop (runId, iterId));"

        "create table if not exists sample_wd_mass"
        "( runId integer not null"
        ", referencedRunId integer not null"
        ", iterId integer not null"
        ", starId text not null"
        ", mass double not null"
        ", clusterMembership double not null"
        ", precursorLogAge double not null"
        ", coolingAge double" // Can sometimes be NaN, which registers as null
        ", logTeff double not null"
        ", logG double not null"
        ", primary key (runId, referencedRunId, iterId, starId)"
        ", foreign key (referencedRunId, iterId) references iteration (runId, id)"
        ", foreign key (referencedRunId, starId) references star (runId, starId)"
        ", foreign key (referencedRunId, iterId) references single_pop (runId, iterId));"
    );
}


template <typename T>
void SqlBackingStore<T>::generateRun()
{
   execOnly(
       "insert into run (id) values (NULL);",
       "Creating run"
       );

    this->run = sqlite3_last_insert_rowid(db.get());
}


template <typename T>
void SqlBackingStore<T>::buildInsertStatement()
{
    dbPrepare("insert into iteration values (?, ?);",
              &insertIteration, "Preparing iteration insert");
}


template <typename T>
Iteration SqlBackingStore<T>::nextIteration()
{
    auto iter = BackingStore<T>::nextIteration();

    sqlite3_bind_int(insertIteration, 1, run);
    sqlite3_bind_int(insertIteration, 2, iter.val);

    dbStepAndReset(insertIteration, "Inserting new iteration");

    return iter;
}


template <typename T>
void SqlBackingStore<T>::execOnly(const string s, const string msg)
{
    dbErrorIf(execOnlyRet(s), msg);
}


template <typename T>
int SqlBackingStore<T>::execOnlyRet(const string s)
{
    int ret;

    do
    {
        ret = sqlite3_exec( db.get(), s.c_str(),
                            nullptr, nullptr, nullptr);

        if ( ret == SQLITE_BUSY )
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    } while ( ret == SQLITE_BUSY );

    return ret;
}


template <typename T>
void SqlBackingStore<T>::dbErrorIf(int ret, const std::string description)
{
    if (ret)
    {
        throw std::runtime_error("\n" + description + " returned: " + std::to_string(ret) + ".\nSQLite says: " + sqlite3_errmsg(db.get()));
    }
}


template <typename T>
void SqlBackingStore<T>::dbStepAndReset(sqlite3_stmt *stmt, const std::string s)
{
    int ret;

    do
    {
        ret = sqlite3_step(stmt);

        if ( ret == SQLITE_BUSY )
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    } while (ret == SQLITE_BUSY );

    if (ret != SQLITE_DONE)
    {
        dbErrorIf(ret, s.c_str());
    }

    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);
}


template <typename T>
void SqlBackingStore<T>::dbPrepare(const std::string sql, sqlite3_stmt **stmt, const std::string msg)
{
    int ret;

    do
    {
        ret = sqlite3_prepare_v2( db.get(), sql.c_str(), -1, stmt, nullptr);

        if ( ret == SQLITE_BUSY )
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    } while ( ret == SQLITE_BUSY );

    dbErrorIf(ret, msg);
}


template <typename T>
void SqlBackingStore<T>::beginTransaction(const string msg)
{
    execOnly(
        "BEGIN TRANSACTION;",
        "Beginning transaction for " + msg
        );
}


template <typename T>
void SqlBackingStore<T>::endTransaction(const string msg)
{
    execOnly(
        "END TRANSACTION;",
        "Ending transaction for " + msg
        );
}


template <typename T>
void SqlBackingStore<T>::execOnlyInTransaction(const string msg, const string sql)
{
    beginTransaction(msg);

    execOnly(sql, msg);

    endTransaction(msg);
}


template <typename T>
void SqlBackingStore<T>::enableForeignKeys()
{
    execOnly(
        "PRAGMA foreign_keys=ON;",
        "Enabling Foreign Keys"
        );
}
