#!/usr/bin/perl

my %tot;
while(my $line=<STDIN>) {
   my @vals=split(/\t/,$line);
   $vals[3] =~ s/^\s+//g;
   my @hits=split(/\s+/,$vals[3]);
   for(my $it = 0; $it < $#hits; $it+=2 ) {
      my $pcnt= $hits[$it+1];
      my $id= $hits[$it];
      if( !$tot{$id} ) {
         $tot{$id}=$pcnt;
      } else {
         $tot{$id} .= ":$pcnt"; 
      }
      #print "debug [$id] $tot{$id}\n";
      if( !$id ) {
         print "huh: $line\n";   
         print "huh1: [$vals[3]]\n";   
         exit(0);
      }
   }
}

foreach my $val (keys %tot) {
   my @vals=split(/:/,$tot{$val});
   my @sval=sort { $a cmp $b } @vals;
   print "$val @sval\n";
}
