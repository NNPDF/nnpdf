#!/bin/bash
# This script generates the binary database from the dump
cd $(dirname $0)
if [ ! -f theory.db ]; then
    sqlite3 theory.db < theory.dat
fi

