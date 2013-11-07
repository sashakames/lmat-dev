#!/usr/bin/perl -w

use strict;

my $startid=10000000;
my %cnt;

while(my $line = <STDIN>) {
   my @vals=split(/\t/,$line);
   $cnt{$vals[0]} .= "$line";
}

my %humanIds;
foreach my $val (keys %cnt) {
   my @vals=split(/\n/,$cnt{$val});
   if( $#vals > 0 ) {
      foreach my $chk (@vals) {
         my $tid=(split(/\t/,$chk))[0];
         my $nid=$tid;
         if( $chk =~ /plasmid/i || $tid == 9606 ) {
            if( $tid == 9606 ) {
               my $hid="";
               if ($chk =~ /chr(\d+|Un|M|X|Y)/ ) {
                  $hid="chr$1";
               }
               if( !$hid ) {
                  #print "error: $chk\n";
                  #exit(0);
                  $nid=9606;
               } elsif( !$humanIds{$hid} ) {
                  $startid++;
                  $nid=$startid;
                  $humanIds{$hid} = $nid;
               } else {
                  $nid = $humanIds{$hid};
               }
            } else {
               $startid++;
               $nid=$startid;
            }
         }
         print "$nid\t$chk\n";
      }
   } else {
         my $tid=(split(/\t/,$vals[0]))[0];
         print "$tid\t$vals[0]\n";
   }
}
