#!/bin/sh

output=$1
#ls -1 ../runtime_inputs/nul*.m9.*gz | perl -ne 'chomp; if(/random.(\d+)./) { $t=$1-19; print "$t ../runtime_inputs/$_\n";}' > cutoffs.kcnt.m9.100.flst
ls -1 ../runtime_inputs/nul*.m9.*gz | perl -ne 'chomp; if(/random.(\d+)./) { $t=$1-19; print "$t $_\n";}' > $output

