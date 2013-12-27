pwd

pushd $1

source $HOME/.bashrc

count=0


for n in `cat complete-sets ` ; do  

    sum=0    
    for val in `cat illum/$n.grp.szs`; do 
	sum=$(( $sum + $val )); 
	
    done

#    if [ $sum -lt 500000000 ] ; then

    echo Found $n with $sum data
    
    subs=`echo $n | awk '{print substr($1, 0, 5);}'`
    
    odir=$1/../work/$subs/$n
    
    echo Dir = $odir
    
    if [ ! -d $odir ] ; then
	mkdir $odir

    
	lst=`cat illum/$n.grp`

	
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

	    popd

	    pwd

	    cmd="./run_lmat.sh --query_file=$path --db_file=/local/ramfs/m9.20mer.16bit.g1000.db --prune_thresh=75 --nullm=no --cs_off --gl_off --threads=40 --odir=$odir --qual_filter=30"
	    echo $cmd
	    sh  $cmd
	   

	    
	done
	time sh do_parallel.sh $odir
	count=$(($count + 1))
	if [ count -eq $2 ] ; then
	    break
	fi
    fi
#    fi

done

