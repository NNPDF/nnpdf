#!/usr/bin/env python

import sqlite3 as lite
import sys,os

# Attempt to find tablulate
import imp
try:
    imp.find_module('tabulate')
    found = True
except ImportError:
    found = False

# Install/import tabulate
if found == False:
    os.system("pip install tabulate --user")
from tabulate import tabulate

# sqlite con
con = None

try:
    con = lite.connect('theory.db')

    cur = con.cursor()
    cur.execute('SELECT SQLITE_VERSION()')

    data = cur.fetchone()

    print("SQLite version: %s" % data)

    cur.execute('SELECT * FROM TheoryIndex')
    col_names = [cn[0] for cn in cur.description]
    col_sub = [col_names[0], col_names[33]]


    table = []
    rows = cur.fetchall()
    for row in rows:
        table.append([row[0], row[36]])

    print(tabulate(table, headers=col_sub))

except lite.Error as e:

    print("Error %s:" % e.args[0])
    sys.exit(1)

finally:

    if con:
        con.close()
