#!/usr/bin/perl -w

use strict;

## simple frequency estimator 

my %save;
my $mrank=$ARGV[0];
my %once;
my %sid;
while(my $line = <STDIN>) {
   chomp($line);
   my ($sum,$id,$pcnt)=(0,"",0);
   if( $line =~ /\t(root|LowScore);(\d+)/ ) {
      ($id,$sum)=($1,$2); 
      $once{$id} = $sum;
   } else {
      my ($ktid,$lineage)=split(/\t/,$line);
      my @ranks=split(/:/,$lineage);
      my $lin="";
      my $cnt=0;
      for(my $it = 0; $it <= $#ranks; $it++) {
         my $rank=$ranks[$it];
         if( $rank =~ /(.*),(.*);(\d+)/ ) {
            my ($rval,$name,$count) = ($1,$2,$3);
            $cnt += $count;
            if( $rval eq $mrank ) {
               $id=$name;
               if(!$once{$id}) {
                  $once{$id} = $cnt;
               } else {
                  $once{$id} += ($cnt-$count);
               }
               last
            }
         }  
      }
   }
}

my $tot=0;

foreach my $lin (keys %once) {
   $tot+= $once{$lin};
}
foreach my $lin (keys %once) {
   my $sum = $once{$lin};
   my $pcnt = $sum / $tot; 
   print "$pcnt $sum $lin\n"; 
}


