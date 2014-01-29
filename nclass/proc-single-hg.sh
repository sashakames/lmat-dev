source $HOME/.bashrc


base_dir=/p/lscratche/ames4/1000humangenome/hg-files
#threads=$4
#dbsrc=$5
subdir=work
#destpath=$6
procset=$1

echo args:  $*



threads=38

dbsrc=/p/lscratche/ames4/m9.20mer.16bit.g1000.db
dbdest=/local/ramfs/m9.20mer.16bit.g1000.db

dbdir=`dirname $dbdest`


if [ ! -f  $dbdest ] ; then
    time cp $dbsrc $dbdest 
fi

hostname

date

subs=`echo $procset | awk '{print substr($base_dir, 0, 5);}'`

subgroup=$base_dir/../$subdir/$subs/

if [ ! -d $subgroup ] ; then
    mkdir $subgroup
fi

    
    
odir=$base_dir/../$subdir/$subs/$procset

rm $odir/*.fa $odir/*.human_matches


echo Dir = $odir

lst=`cat $base_dir/illum/$procset.grp`

	
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
	   
	    time sh do_parallel.sh $odir $path

	done

inf=$odir/$procset.new.hg
outf=$odir/$procset.jf
xcover=30

nthreads=$(( threads + 2 ))

time $HOME/bin/jellyfish count -m 20 -o $outf -t $nthreads -c 6 -s 100000000 -C $inf 
$HOME/bin/jellyfish dump -c -L $xcover $outf | awk '{print $1;}' | sort | uniq > $outf.$xcover.kmerdump

echo $1 has completed

date