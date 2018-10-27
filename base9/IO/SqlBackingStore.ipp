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
        #include "sql/BuildTables.sql"
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
    return doWhileBusy([this, &s]
                       {
                           return sqlite3_exec( db.get(), s.c_str(),
                                                nullptr, nullptr, nullptr);
                       });
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
int SqlBackingStore<T>::doWhileBusy(std::function<int()> func)
{
    int ret;

    do
    {
        ret = func();

        if ( ret == SQLITE_BUSY )
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    } while (ret == SQLITE_BUSY );

    return ret;
}

template <typename T>
void SqlBackingStore<T>::dbStepAndReset(sqlite3_stmt *stmt, const std::string s)
{
    int ret = doWhileBusy([stmt]() { return sqlite3_step(stmt); });

    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);
}


template <typename T>
void SqlBackingStore<T>::dbPrepare(const std::string sql, sqlite3_stmt **stmt, const std::string msg)
{
    int ret = doWhileBusy([this, &sql, stmt]()
                          {
                              return sqlite3_prepare_v2( db.get(),
                                                         sql.c_str(),
                                                         -1,
                                                         stmt,
                                                         nullptr);
                          });

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
