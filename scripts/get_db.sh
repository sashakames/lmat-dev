if test $# = 0
then
    echo "$0: LMAT database auto-download utility"
    echo ""
    echo ""
    echo "Usage: $0 <database name|inputs name> [Destination path for database/input files]"
    echo "Example databases"
    echo "kML.v4-14.20.g10.db  - Microbial marker database (small database for fast microbial profiling"
    echo "kML+Human.v4-14.20.g10.db  - Microbial marker database with explicit human read tagging (small database for fast microbial profiling)"
    echo "lmat-4-14.20mer.db  - Fullsized database for extensive read binning"
    echo "lmat.genes.7-14.db  - Gene database for gene name binning"
    echo "lmat-world-region  - Database for binning human reads by world region"
    echo ""
    echo "For all databases except world-region, use '$0 inputs 04072014 to download the needed auxilliary data"
    echo "For world-region use '$0 inputs world-region"
    echo ""
    echo "Legacy usage for older databases"
    echo "Usage: $0 <kML-18mer-large|kML-18mer-medium|kML-18mer-small|gene-20mer|kFull-20mer|inputs> [Destination path for database/input files]"
    echo ""
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
echo Download complete.  When running LMAT set --db_file=$outpath/$dbname.db

elif [ $dbname == "kML-18mer-small" ]
then
dbname=kML.18mer.16bit.reduced
echo downloading...
wget -q -O -  ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/$dbname.db.gz | gunzip -c > $outpath/$dbname.db
echo Download complete.  When running LMAT set --db_file=$outpath/$dbname.db

elif [ $dbname == "gene-20mer" ]
then
dbname=gene.20mer
for suffix in a b c d e
do
wget -q -O -  ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/GeneDB/$dbname.db.a$suffix.gz | gunzip -c >> $outpath/$dbname.db
echo "Part $suffix out of 5 done"
done
echo Download complete.  When running LMAT set --db_file=$outpath/$dbname.db
elif [ $dbname == "kFull-20mer" ]
then

   mx=19
   for suffix in `seq 0 $mx` ; do
      file=kFull.20mer.g1000.part.$suffix.lzma
      echo "Retrieve $file"
      wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/20merFullDB/$file | unlzma -c >> $outpath/m9.20mer.16bit.g1000.db
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
echo Download complete.  When running LMAT set --db_file=$outpath/kML.18mer.no_prune.16bit.db 

elif [ $dbname == "inputs" ]
then
echo "Downloading LMAT runtime-input files to $outpath"
input_file=$2
wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/runtime_inputs/$input_file.tgz | tar -C $outpath/ -zxvf - 
abspath=`readlink -f $outpath`
echo "For LMAT to run correctly, please set the LMAT_DIR environment variable to $abspath"

else
   ## Now assume naming convention to avoid updating this file for evry new database
   wget -q -O $outpath/dbinfo ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/$dbname/dbinfo 
   mx=`head -1 $outpath/dbinfo | cut -f1`
   cmprs=`head -1 $outpath/dbinfo | cut -f2`
   mbytes=`head -1 $outpath/dbinfo | cut -f3`
   rm -f $outpath/dbinfo
   echo "Debug: $mx $cmprs $mbytes"
   if [ $mx == -1 ] ; then
      file=$dbname.$cmprs
      if [ $cmprs == "lzma" ] ; then
         wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/$dbname/$file | unlzma > $outpath/$dbname
      elif [ $cmprs == "gz" ] ; then 
         wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/$dbname/$file | gunzip > $outpath/$dbname
      else 
         echo "Unrecognized compression, failed to download"
      fi
   else 
      for suffix in `seq 0 $mx` ; do
         file=$dbname.$suffix.$cmprs
         echo "Retrieve $file"
         if [ $cmprs == "lzma" ] ; then
            wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/$dbname/$file | unlzma -c >> $outpath/$dbname
         elif [ $cmprs == "gz" ] ; then 
            wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/$dbname/$file | gunzip -c >> $outpath/$dbname
         else 
            echo "Unrecognized compression, failed to download"
         fi
         echo part $suffix out of $mx done
         size=`stat $outpath/$dbname | grep Size | awk '{print $2}'`
         if [ $size -gt $mbytes ] ; then
            truncate -s $mbytes $outpath/$dbname
            break
         fi
      done
      echo Download complete.  When running LMAT set --db_file=$outpath/$dbname
   fi

fi

echo "Download complete"
