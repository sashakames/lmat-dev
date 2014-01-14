pwd

source $HOME/.bashrc

for n in `cat $1/complete-sets ` ; do  

    sum=0    
    for val in `cat $1/illum/$n.grp.szs`; do 
	sum=$(( $sum + $val )); 
	
    done

#    if [ $sum -lt 500000000 ] ; then

    echo Found $n with $sum data
    
    subs=`echo $n | awk '{print substr($1, 0, 5);}'`

    subdir=$1/../work/$subs/ 
    
    if [ ! -d $subdir ] ; then
	mkdir $subdir
    fi

    odir=$subdir/$n


    
    echo Dir = $odir
    
    if [ ! -d $odir ] ; then
	mkdir $odir

    
	lst=`cat $1/illum/$n.grp`

	
	for m in $lst  ; do
	    
            fn=`basename $m`

            p1=/p/lscratche/allen99/1000genomes/$fn
            p2=/p/lscratche/allen99/1000genomes/test/$fn
	    
	    if [ -f $p1 ] ; then  
		path=$p1
            else
		path=$p2
            fi
	    
	    echo $path

	    cmd="./run_lmat.sh --query_file=$path --db_file=/local/ramfs/m9.20mer.16bit.g1000.db --prune_thresh=75 --nullm=no --cs_off --gl_off --threads=40 --odir=$odir --qual_filter=30 --min_read_kmer=1"
	    echo $cmd
	    sh  $cmd
	   
	    time sh do_parallel.sh $odir $path

	    
	done



    fi
#    fi

done

