#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "sqlite/sqlite3.h"


class SqlQueryRunner
{
  public:
    SqlQueryRunner(std::string filename, std::string sql)
    {
        sqlite3_open(filename.c_str(), &db);

        while ( SQLITE_BUSY == sqlite3_prepare_v2( db, sql.c_str(), -1,
                                                   &stmt, nullptr ))
        { ; }
    }

    ~SqlQueryRunner()
    {
        sqlite3_finalize(stmt);
        stmt = nullptr;

        sqlite3_close(db);
        db = nullptr;
    }

    std::vector<std::vector<std::string>> run()
    {
        int ret;

        do
        {
            ret = sqlite3_step(stmt);
        } while (ret == SQLITE_BUSY);

        if ( ret == SQLITE_DONE )
        {
            std::cout << "Query returned no values" << std::endl;
        }
        else if (ret == SQLITE_ROW)
        {
            auto columns = sqlite3_column_count(stmt);

            table.emplace_back();

            for (int i = 0; i < columns; ++i)
            {
                table.back().push_back(std::string(sqlite3_column_name(stmt, i)));
                widths.push_back(table.back().back().size());
            }


            while ( ret == SQLITE_ROW )
            {
                table.emplace_back();

                for (int i = 0; i < columns; ++i)
                {
                    std::string s = std::string(reinterpret_cast<const char*>
                                                    (sqlite3_column_text(stmt, i)));
                    table.back().push_back(s);

                    auto size = s.size();

                    if (size > widths.at(i))
                    {
                        widths.at(i) = size;
                    }
                }

                ret = sqlite3_step(stmt);
            }
        }
        else
        {
            std::cerr << "Query returned code: " << std::to_string(ret) << ".\nSQLite says: " << sqlite3_errmsg(db) << std::endl;
        }

        sqlite3_reset(stmt);

        return table;
    }


    void prettyPrint()
    {
        for (size_t i = 0; i < table.size(); ++i)
        {
            const auto &v = table.at(i);

            if ( i == 1 )
            {
                for (size_t j = 0; j < v.size(); ++j)
                {
                    std::cout << std::string(widths.at(j), '-') << "  ";
                }

                std::cout << std::endl;
            }

            for (size_t j = 0; j < v.size(); ++j)
            {
                if ( j > 0 )
                {
                    std::cout << "  ";
                }

                std::cout << std::setw(widths.at(j)) << std::left << v.at(j);
            }

            std::cout << std::endl;
        }
    }

  private:
    sqlite3 *db = nullptr;
    sqlite3_stmt *stmt = nullptr;

    std::vector<std::vector<std::string>> table;

    std::vector<size_t> widths;
};
