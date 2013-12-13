if test $# = 0
then
    echo "$0: LMAT database auto-download utility"
    echo "Usage: $0 <kML-18mer-large|kML-18mer-medium|kML-18mer-small|gene-20mer|kFull-20mer|inputs> [Destination path for database/input files]"

    echo
    echo "'inputs' downloads and extracts the LMAT 'runtime_inputs' files."
    echo "Please see LMAT documentation for more details."
    exit
fi

if test $# = 1
then
outpath="."
else
outpath=$2
fi

dbname=$1

if ! test -d $outpath; then
    echo "$outpath not found, creating directory"
    mkdir $outpath
    if ! test -d $outpath; then
	echo Could not create directory for database download.  Please ensure that the parent directory name is correct.
	exit 1
    fi
fi


if [ $dbname == "kML-18mer-medium" ]
then
dbname=kML.18mer.16bit
for suffix in a b c d e
do
    wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/$dbname.db.$suffix.gz | gunzip -c >> $outpath/$dbname.db
echo "Part $suffix out of 5 done"
done
echo Download complete.  When running LMAT set --marker_library=$outpath/$dbname.db

elif [ $dbname == "kML-18mer-small" ]
then
dbname=kML.18mer.16bit.reduced
echo downloading...
wget -q -O -  ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/$dbname.db.gz | gunzip -c > $outpath/$dbname.db
echo Download complete.  When running LMAT set --marker_library=$outpath/$dbname.db

elif [ $dbname == "gene-20mer" ]
then
dbname=gene.20mer
for suffix in a b c d e
do
wget -q -O -  ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/GeneDB/$dbname.db.a$suffix.gz | gunzip -c >> $outpath/$dbname.db
echo "Part $suffix out of 5 done"
done
echo Download complete.  When running LMAT set --genedb_file=$outpath/$dbname.db

elif [ $dbname == "kFull-20mer" ]
then

   mx=19
   for suffix in `seq 0 $mx` ; do
      file=kFull.20mer.g1000.part.$suffix.lzma
      echo "Retrieve $file"
      wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/20merFullDB/$file | unlzma -c >> $outpath/m9.20mer.16bit.g1000.db
      #echo part $(( 1 + $suffix )) out of $mx done
      echo part $suffix out of $mx done
      size=`stat $outpath/m9.20mer.16bit.g1000.db | grep Size | awk '{print $2}'`
      if [ $size -gt 400000000000 ] ; then
         truncate -s 400GB $outpath/m9.20mer.16bit.g1000.db
         break
      fi
   done
   echo Download complete.  When running LMAT set --db_file=$outpath/m9.20mer.16bit.g1000.db 

elif [ $dbname == "kML-18mer-large" ]
then
for suffix in `seq 0 7` ; do
wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/kML.18mer.no_prune.16bit.part.$suffix.lzma | unlzma -c >> $outpath/kML.18mer.no_prune.16bit.db
echo part $(( 1 + $suffix )) out of 8 done
done
echo Download complete.  When running LMAT set --marker_libary=$outpath/kML.18mer.no_prune.16bit.db 

elif [ $dbname == "inputs" ]
then
echo "Downloading LMAT runtime-input files to $outpath"
wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/runtime_inputs/09042013.tar.gz | tar -C $outpath/ -zxvf - 
abspath=`readlink -f $outpath`
echo "For LMAT to run correctly, please set the LMAT_DIR environment variable to $abspath"

else
echo "Unrecognized option: [$dbname]" 
exit 0
fi

echo "Download complete"
