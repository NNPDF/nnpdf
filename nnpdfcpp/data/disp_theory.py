#!/usr/bin/python

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
    
    print "SQLite version: %s" % data    

    cur.execute('SELECT * FROM TheoryIndex')
    col_names = [cn[0] for cn in cur.description]

    table = []

    rows = cur.fetchall()

    print tabulate(rows, headers=col_names)

    theoryID = int(raw_input("Please select a table ID: "))
    if 0 <= theoryID <= len(rows)-1:
        os.system("wget http://pcteserver.mi.infn.it/~apfelcomb/commondatatheory/theory_%d.tgz" % theoryID)
        os.system("tar -xvzf theory_%d.tgz" % theoryID)
        os.system("rm theory_%d.tgz" % theoryID)
    else:
        print "Invalid theory ID: " + str(theoryID)
        sys.exit(1)
    
except lite.Error, e:
    
    print "Error %s:" % e.args[0]
    sys.exit(1)
    
finally:
    
    if con:
        con.close()
