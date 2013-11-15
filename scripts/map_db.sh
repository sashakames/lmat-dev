if test $# = 0
then
    echo "$0: LMAT database auto-download with memory mapping utility"
    echo "Usage: $0 <kML.18mer.full> [Destination path for database/input files]"
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


if [ $dbname == "kML.18.full" ]
then
for suffix in 0 1 2 3 4 5 6 7
do
    wget -q -O - ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/18merML/$dbname.tax_histo.$suffix.gz | gunzip -c > $outpath/$dbname.$suffix
    echo "Part $suffix out of 7 done"
    echo "$outpath/$dbname.db" >> $outpath/$dbname.build_lst.txt
done
$LMAT_DIR/../bin/make_db_table -f 16 -i $outpath/$dbname.build_lst.txt -o $outpath/$dbname.db -k 18 -s 75 -h 2000000000

else
echo "Unrecognized option: [$dbname]" 
exit 0
fi

echo "Download complete"
