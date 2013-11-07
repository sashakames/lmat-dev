CMD=./read_search
DBFILE=../src/kmerdb/examples/tests/data/test.fa.int.bin.no_loc.12
READFILE=../src/kmerdb/examples/tests/data/reads.fa
OUTDIR=../src/kmerdb/examples/tests/data/
NTHREADS=4
WHOLECMD="$CMD -d $DBFILE -i $READFILE -t $NTHREADS -o $OUTDIR"
echo $WHOLECMD
$WHOLECMD

