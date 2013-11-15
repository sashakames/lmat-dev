#!/bin/sh

lpe=2 #lines_per_entry 2 for fasta and 4 for fastq

PATH_TO_INPUT_DATA=./orig
for f in $PATH_TO_INPUT_DATA/*fasta; do
    n=${f%.fasta}
    n=${n##*/}

#    fN=`wc -l $f|sed 's/ .*//'`
    echo $n
    echo $f
    for N in 1000 9000 27000 81000 2187000 ; do

        let "Q=$N*$lpe"
        echo "N=$N Q=$Q fN=$fN"
	
	headin=$n.$N.head
	tailin=$n.$N.tail

	if [ ! -f $headin ] ; then
	    echo head $Q $f
	    head -n $Q  $f > $headin
	    
	fi

	if [ ! -f $tailin ] ; then
	    echo tail $Q $f
	    tail -n $Q $f > $tailin	    
	fi


	for thresh in 20 30 40 50 75 100; do

            odir=dir.$thresh.$n.$1.$N/
	    
            if [ ! -d $odir ] ; then
		echo creating $odir
		mkdir $odir
            fi
	    echo $odir

	    
#            run_lmat.sh --query_file=$headin --marker_library=/local/ramfs/$1 --marker_min_score=0 --gl_off --odir=$odir --min_read_kmer=5 --cs_off
#            run_lmat.sh --query_file=$tailin --marker_library=/local/ramfs/$1 --marker_min_score=0 --gl_off --odir=$odir --min_read_kmer=5 --cs_off
	    

        done
        
    done

done
