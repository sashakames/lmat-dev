#!/usr/bin/perl -w

use strict;

my %save;
my $file=$ARGV[0];
open(FILE,$file)|| die "fail $file\n";
 
while(my $line = <FILE>) {
   chomp($line);
   my @vals=split(/\s+/,$line);
   #print "save [$vals[0]]\n";
   $save{$vals[0]}=1;
}
close(FILE);
while(my $line = <STDIN>) {
   my @vals=split(/\s+/,$line);
   if( $save{$vals[0]} ) {
      print $line;
   }
}
