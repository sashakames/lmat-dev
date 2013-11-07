#!/bin/sh -xvf

ddir=/usr/mic/post1/metagenomics
tax_histo -d $ddir/small_test_data/all_virus.int.bin.no_loc.16 -t $ddir/taxonomy/tax_nodes_ncbi.dat -i $ddir/taxonomy/viral_gid2tid.dat
