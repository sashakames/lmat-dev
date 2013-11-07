if test $# = 0
then
    echo "$0: LMAT database auto-download utility"
    echo "Usage: $0 <kML.18mer.16bit|kML.18mer.16bit.reduced|m9.0904.16|gene.20mer|inputs> [Destination path for database/input files]"
    echo
    echo "'inputs' downloads and extracts the LMAT 'runtime_inputs' files."
    echo "Please see LMAT documentation for more details."
    exit
fi

if test $# = 1
then
outpath="."
# assume no u
else
outpath=$2
fi

dbname=$1


if [ $dbname == "kML.18mer.16bit" ]
then
for suffix in a b c d e
do
    wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/$dbname.db.$suffix.gz | gunzip -c >> $outpath/$dbname.db
echo "Part $suffix out of 5 done"
done
elif [ $dbname == "kML.18mer.16bit.reduced" ]
then
echo downloading...
wget -q -O -  ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/$dbname.db.gz | gunzip -c > $outpath/$dbname.db

elif [ $dbname == "gene.20mer" ]
then
for suffix in a b c d e
do
wget -q -O -  ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/GeneDB/$dbname.db.a$suffix.gz | gunzip -c >> $outpath/$dbname.db
echo "Part $suffix out of 5 done"
done
elif [ $dbname == "m9.0904.16" ]
then
echo "Database coming soon!"

elif [ $dbname == "inputs" ]
then
echo "Downloading LMAT runtime-input files to $outpath"
#wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/runtime_inputs/09042013.tar.gz | tar -C $outpath/ -zxvf - 
abspath=`readlink -f $outpath`
echo "For LMAT to run correctly, please set the LMAT_DIR environment variable to $abspath"

else
echo "Unrecognized option: [$dbname]" 
exit 0
fi

echo "Download complete"
