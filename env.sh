#!/bin/sh
export LD_LIBRARY_PATH=$PJPATH/lib:$LD_LIBRARY_PATH
export METAG_DIR=`pwd`/..
export PJPATH=$METAG_DIR/third-party/perm-je-0.9.3

export WITH_PJMALLOC=1
#export WITH_PJMALLOC=0
export BOOST=0
export USE_SORTED_DB=1

