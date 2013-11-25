#!/bin/sh

export LMAT_DIR=/usr/mic/post1/metagenomics/runtime_inputs/09042013/

export METAG_DIR=$HOME/LMAT-1.2.1/

bin_dir=$METAG_DIR/bin

export PATH=$PATH:$bin_dir




lst=`cat $1`

cd $bin_dir

for fn in $lst ; do

    dir=`dirname $fn`
    tmp=`basename $fn`
    n=${tmp%.*}
    
    
   
#for thresh in 20 30 40 50 75 100; do


    # use thresh
    odir=$dir/dir.$n.$2/
	    
    if [ ! -d $odir ] ; then
	echo creating $odir
	mkdir $odir
    fi
    echo $odir

    
  
    sh run_lmat.sh --query_file=$fn --db_file=/local/ramfs/m9.20mer.16bit.g1000.db --gl_off --odir=$odir --cs_off --threads=80 --min_read_kmer=5
#            run_lmat.sh --query_file=$tailin --marker_library=/local/ramfs/$1 --marker_min_score=0 --gl_off --odir=$odir --min_read_kmer=5 --cs_off
	    
#done

done
