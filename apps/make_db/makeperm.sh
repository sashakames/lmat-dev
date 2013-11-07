#!/bin/sh
make WITH_PJMALLOC=1 PJPATH=/home/gokhale2/local make_db
make WITH_PJMALLOC=1 PJPATH=/home/gokhale2/local read_db

