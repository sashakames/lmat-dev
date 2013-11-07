export LD_LIBRARY_PATH=/usr/mic/post1/metagenomics/lib
./kmerPrefixCounterSort -i /p/lscratche/hysom/m9.09042013.fa.int -o  /p/lscratche/ames4/m9.09042013-kmers/m9.09042013.$1 -k 20 -l 3 -f $1
