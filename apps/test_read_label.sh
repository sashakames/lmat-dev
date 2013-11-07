#!/bin/sh -xvf

#remmember /usr/mic/post1/ is nsf - so this is for debuggin
ddir=/usr/mic/post1/metagenomics
qfile=/usr/mic/post1/metagenomics/metagenomes/SRX022172.fasta
metag_dir=/g/g21/allen99/repo/metag_repo
bin=$metag_dir/dev/apps
tdir=$metag_dir/taxonomy/kpath

read_label -t 1 -i $qfile -d $ddir/tmp_kmer_dbs/all_virus.int.bin.no_loc.16.flst -c $tdir/kpath_taxonomy.dat -o debug.txt
