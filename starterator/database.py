#!/usr/bin/env python
#
# Copyright 2009 Facebook
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may obtain
# a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.


try:
    import pymysql
    pymysql.install_as_MySQLdb()
    import pymysql as MySQLdb
except ImportError:
    pymysql = None
    MySQLdb = None
import sqlite3
import time
import os
from .utils import get_config, StarteratorError


class _MySQLOperationalErrorPlaceholder(Exception):
    pass


MYSQL_OPERATIONAL_ERROR = (
    MySQLdb.OperationalError
    if MySQLdb is not None
    else _MySQLOperationalErrorPlaceholder
)

_PROCESS_DB = None
_PROCESS_DB_PID = None

class DB(object):

    def __init__(self):
        self._last_use_time = time.time()
        self.max_idle_time = float(1800)
        self._db = None
        self.backend = "mysql"
        self.target = ""

        sqlite_path = os.getenv("DB_SQLITE_PATH", "").strip()
        if sqlite_path:
            self.backend = "sqlite"
            self._sqlite_path = os.path.abspath(os.path.expanduser(sqlite_path))
            self.target = self._sqlite_path
            self.host = self._sqlite_path
            self._db_args = {}
        else:
            if MySQLdb is None:
                raise StarteratorError(
                    "PyMySQL is required for MySQL mode. Install pymysql or set DB_SQLITE_PATH for SQLite mode."
                )
            args = {}
            config = get_config()
            args["user"] = os.getenv("DB_USER", config.get("database_user", ""))
            args["password"] = os.getenv("DB_PASSWORD", config.get("database_password", ""))
            args["database"] = os.getenv("DB_NAME", config.get("database_name", ""))
            args["host"] = os.getenv("DB_HOST", config.get("database_server", "localhost"))
            args["port"] = int(os.getenv("DB_PORT", config.get("database_port", "3306")))

            # Only use unix_socket for local connections
            if args["host"] in ["localhost", "127.0.0.1", "::1"] and not os.getenv("DB_HOST"):
                args["unix_socket"] = "/var/run/mysqld/mysqld.sock"
            self.host = "{0}:{1}".format(args['host'], args['port'])
            self.target = self.host
            self._db_args = args

        try:
            self.reconnect()
        except Exception:
            raise

    def __del__(self):
        self.close()

    def close(self):
        """Closes this database connection."""
        if getattr(self, "_db", None) is not None:
            self._db.close()
            self._db = None

    def reconnect(self):
        """Closes the existing database connection and re-opens it."""
        self.close()
        if self.backend == "sqlite":
            if not os.path.exists(self._sqlite_path):
                raise StarteratorError("SQLite database file not found: %s" % self._sqlite_path)
            try:
                self._db = sqlite3.connect(self._sqlite_path)
            except sqlite3.Error as e:
                raise StarteratorError("Error connecting to SQLite database at %s: %s" % (self._sqlite_path, e))
        else:
            self._db = MySQLdb.connect(**self._db_args)

    def iter(self, query, params):
        """Returns an iterator for the given query and parameters."""
        self._ensure_connected()
        if self.backend == "mysql":
            cursor = MySQLdb.cursors.SSCursor(self._db)
        else:
            cursor = self._db.cursor()
        try:
            self._execute(cursor, query, params)
            for row in cursor:
                yield row
        finally:
            cursor.close()

    def query(self, query, params=None):
        '''
        cursor = self._cursor()
        try:
            self._execute(cursor, query, params)
            result = cursor.fetchall()
            return result.
        except:
            self.reconnect()
            self.query(query, params)
        finally:
            cursor.close()
        '''
        #Potential fix for connection error from github copilot
        cursor = self._cursor()
        try:
            self._execute(cursor, query, params)
            result = cursor.fetchall()
            return result
        except MYSQL_OPERATIONAL_ERROR as e:
            if self.backend != "mysql":
                raise
            print("OperationalError: %s" % e)
            self.reconnect()
            return self.query(query, params)
        finally:
            if cursor:
                try:
                    cursor.close()
                except Exception:
                    pass


    def get(self, query, params):
        """Returns the first row returned for the given query."""
        rows = self.query(query, params)
        if not rows:
            return None
        elif len(rows) > 1:
            raise Exception("Multiple rows returned for Database.get() query")
        else:
            return rows[0]

    def _ensure_connected(self):
        # Mysql by default closes client connections that are idle for
        # 8 hours, but the client library does not report this fact until
        # you try to perform a query and it fails.  Protect against this
        # case by preemptively closing and reopening the connection
        # if it has been idle for too long (7 hours by default).
        if self.backend == "sqlite":
            if self._db is None:
                self.reconnect()
            self._last_use_time = time.time()
            return
        if (self._db is None or (time.time() - self._last_use_time > self.max_idle_time)):
            self.reconnect()
        self._last_use_time = time.time()

    def _cursor(self):
        self._ensure_connected()
        return self._db.cursor()

    def _execute(self, cursor, query, params):
        if self.backend == "sqlite":
            sqlite_query, sqlite_params = self._convert_query_for_sqlite(query, params)
            try:
                if sqlite_params:
                    return cursor.execute(sqlite_query, sqlite_params)
                return cursor.execute(sqlite_query)
            except sqlite3.Error as e:
                raise StarteratorError(
                    "SQLite query failed on %s: %s. Query: %s" % (self.target, e, query)
                )

        try:
            if params is None:
                return cursor.execute(query)
            elif isinstance(params, tuple):
                return cursor.execute(query, params)
            else:
                return cursor.execute(query, (params,))
        except MYSQL_OPERATIONAL_ERROR:
            print("Error connecting to MySQL on %s" % self.host)
            self.close()
            raise StarteratorError("Error connecting to database! Please enter correct login credentials in Preferences menu.")

    def _normalize_params(self, params):
        if params is None:
            return []
        if isinstance(params, (tuple, list)):
            return list(params)
        return [params]

    def _convert_query_for_sqlite(self, query, params):
        param_values = self._normalize_params(params)
        converted_query_parts = []
        sqlite_params = []
        param_index = 0
        i = 0

        while i < len(query):
            char = query[i]
            if char == '%' and i + 1 < len(query):
                marker = query[i + 1]
                if marker == '%':
                    converted_query_parts.append('%')
                    i += 2
                    continue
                if marker == 's':
                    if param_index >= len(param_values):
                        raise StarteratorError(
                            "SQLite parameter mismatch: query expects more parameters. Query: %s" % query
                        )
                    value = param_values[param_index]
                    param_index += 1
                    if isinstance(value, (tuple, list)):
                        if len(value) == 0:
                            raise StarteratorError(
                                "SQLite does not support empty parameter lists for IN clauses. Query: %s" % query
                            )
                        converted_query_parts.append("(" + ",".join(["?"] * len(value)) + ")")
                        sqlite_params.extend(value)
                    else:
                        converted_query_parts.append("?")
                        sqlite_params.append(value)
                    i += 2
                    continue
            converted_query_parts.append(char)
            i += 1

        if param_index != len(param_values):
            raise StarteratorError(
                "SQLite parameter mismatch: %d unused parameter(s). Query: %s" % (
                    len(param_values) - param_index, query
                )
            )

        return "".join(converted_query_parts), tuple(sqlite_params)
            
class Row(dict):
    """A dict that allows for object-like property access syntax."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)


def get_db():
    global _PROCESS_DB, _PROCESS_DB_PID
    current_pid = os.getpid()

    if _PROCESS_DB is None or _PROCESS_DB_PID != current_pid:
        if _PROCESS_DB is not None:
            try:
                _PROCESS_DB.close()
            except Exception:
                pass
        _PROCESS_DB = DB()
        _PROCESS_DB_PID = current_pid

    return _PROCESS_DB
