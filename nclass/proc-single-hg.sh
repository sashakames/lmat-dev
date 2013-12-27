source $HOME/.bashrc

subs=`echo $2 | awk '{print substr($1, 0, 5);}'`
    
odir=$1/../work/$subs/$2

echo Dir = $odir

lst=`cat $1/illum/$2.grp`

	
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

	    pwd

	    cmd="./run_lmat.sh --query_file=$path --db_file=/local/ramfs/m9.20mer.16bit.g1000.db --prune_thresh=75 --nullm=no --cs_off --gl_off --threads=40 --odir=$odir --qual_filter=30"
	    echo $cmd
	    sh  $cmd
	   

	    
	done
	time sh do_parallel.sh $odir