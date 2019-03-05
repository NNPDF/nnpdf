#!/bin/bash
# This script dumps the theory database into a text file
sqlite3 theory.db ".dump" > theory.dat
