#!/usr/bin/perl
#
use strict;

my @private_ids=("Calvin Keeler", "CDC", "Eureka", "Imigene", "Keim NAU", "LANL", "LLNL", "NBACC", "NBFAC", "PIADC", "USAMRIID", "USDA", "UTMB", "UVIC");

while(my $line =<STDIN>) {
   my @vals=split(/\t/,$line);
   ## sequences without GenBank ID need
   ## not be proprietary if they're assemblies/concatenations
   ## of draft contigs
   my $private=0;
   if( $vals[0] >= 10000000) {
      ## no ncbi taxonomy was assigned so skip
      $private=1;
   }
   elsif( $vals[1] == -1 ) {
      foreach my $id (@private_ids) {
         if( $line =~ / $id / ) {
            $private=1;
            last;
         }
      }
   }
   if( !$private ) {
      print $line;
   }
}
