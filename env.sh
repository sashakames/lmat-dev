#!/bin/sh
export METAG_DIR=`pwd`
export PJPATH=/usr/mic/post1/metagenomics/perm-je-0.9.7
export LD_LIBRARY_PATH=$PJPATH/lib:$LD_LIBRARY_PATH

export WITH_PJMALLOC=1
#export WITH_PJMALLOC=0
export BOOST=0
export USE_SORTED_DB=1

export TAXID_SIZE=16
