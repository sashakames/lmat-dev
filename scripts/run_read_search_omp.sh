#!/bin/sh

LVALS="4 8 12 16 20 24 32 40 64 80 120"

TRIALS=4

CMD=$HOME/metag_repo/dev/apps/read_search

OUTDIR=/p/lscratchc/ames4/rso-1

OUTBASE=$OUTDIR/out
TIMEBASE=$OUTDIR/time

DBFILE=/usr/mic/post1/metagenomics/small_test_data/all_virus.int.bin.no_loc.18
READFILE=/usr/mic/post1/metagenomics/metagenomes/SRX022172.fasta


echo "command: $CMD" > $OUTDIR/argsfile.$HOSTNAME
echo "thread_vals: $LVALS" >> $OUTDIR/argsfile.$HOSTNAME
echo "trials: $TRIALS" >> $OUTDIR/argsfile.$HOSTNAME
echo "dbfile: $DBFILE" >> $OUTDIR/argsfile.$HOSTNAME
echo "readfile: $READFILE" >> $OUTDIR/argsfile.$HOSTNAME

chown $LOGNAME $OUTDIR/argsfile

for i in $LVALS 
 do
 for (( j = 1; j <= $TRIALS; j++ ))
   do
   echo 3 > /proc/sys/vm/drop_caches
#   time ($CMD -d $DBFILE -i $READFILE -t $i -o /tmp/rso > $OUTBASE.$i.$j.$HOSTNAME ) 2> $TIMEBASE.$i.$j.$HOSTNAME
   time ($CMD -d $DBFILE -i $READFILE -t $i -o /tmp/rso -r 1 > $OUTBASE.$i.$j.$HOSTNAME ) 2> $TIMEBASE.$i.$j.$HOSTNAME
   chown ames4 $TIMEBASE.$i.$j.$HOSTNAME
   chown ames4 $OUTBASE.$i.$j.$HOSTNAME
done

done


