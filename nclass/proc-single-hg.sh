source $HOME/.bashrc


base_dir=$1
threads=$4
dbsrc=$5
subdir=$2
destpath=$6
procset=$3

echo args:  $*


if [ ! -m $threads ] ; then
    threads=40
fi

dbdest=/local/ramfs/m9.20mer.16bit.g1000.db

# if [ -f $dbsrc ] ; then

#     dbname=`basename $dbsrc`
#     dbdest=$destpath/$dbname
#     source $HOME/.bashrc
    
#     if [ ! -f  $dbdest ] ; then
# 	time cp $dbsrc $dbdest 
#     fi

# else    

# fi 



subs=`echo $procset | awk '{print substr($1, 0, 5);}'`

subgroup=$base_dir/../$subdir/$subs/

if [ ! -d $subgroup ] ; then
    mkdir $subgroup
fi

    

odir=$subgroup/$n
    
echo Dir = $odir
    
if [ ! -d $odir ] ; then
    mkdir $odir
fi

    
odir=$basedir/../$subdir/$subs/$procset

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

	    cmd="./run_lmat.sh --query_file=$path --db_file=$dbdest --prune_thresh=75 --nullm=no --cs_off --gl_off --threads=$threads --odir=$odir --qual_filter=30 --min_read_kmer=1"
	    echo $cmd
	    sh  $cmd
	   

	    
	done
	time sh do_parallel.sh $odir