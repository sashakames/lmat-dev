#!/usr/bin/env perl
use strict;

my %dict;

while(my $fn = <STDIN> ) {
   chomp($fn);
   my @vals=split(/ /,$fn);
   foreach my $val ( @vals ) { 
      if(!$val) {
         last;
      }
      open(FILE,$val) || die "fail $val\n";

      while(my $line=<FILE>) { 
         chomp($line);
         my @parts = split(/;/,$line);
         my $cnt = $parts[$#parts];
         my $id=$parts[0];
         #print "debug: $id $cnt\n";
         if( $dict{$id} ) {
            $dict{$id} += $cnt;
         } else {
            $dict{$id} = $cnt;
         }
      }
      close(FILE);
   }
}
foreach my $k (keys %dict) {
   print "$dict{$k}\t$k\n";
}
