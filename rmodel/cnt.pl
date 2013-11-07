#!/usr/bin/perl
#
my %tot;
while(my $line=<STDIN>) {
   next if( $line !~ /^read_label/);
   my @vals=split(/ /,$line);
   my $cov=$vals[2];
   for(my $it = 3; $it < $#vals - 2; $it+=2 ) {
      my $pcnt= $vals[$it+1];
      my $id= $vals[$it];
      if( !$tot{$id} ) {
         $tot{$id}=$pcnt;
      } else {
         my $t=$pcnt;
         $tot{$id} .= ":$t"; 
      }
   }
}

foreach my $val (keys %tot) {
   my @vals=split(/:/,$tot{$val});
   my @sval=sort { $a cmp $b } @vals;
   print "$val @sval\n";
}
