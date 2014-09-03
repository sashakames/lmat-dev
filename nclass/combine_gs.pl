#!/usr/bin/env perl 

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use strict;

my $rankFile=$ARGV[0];
my %rcnt;
my %score;
my %sid;
while(my $fn = <STDIN>) {
   chomp($fn);
   my @vals=split(/ /,$fn);
   foreach my $val (@vals) {
      if(!$val ) {
         last;
      }
      open(FILE,$val) || die "fail $val\n";
      while(my $line = <FILE> ) {
        chomp($line);
        my @parts = split(/\t/,$line); 
        my $id="$parts[2] $parts[3]";
        #$sid{$id}="$parts[3]\t$parts[4]\t$parts[5]\t$parts[6]\t$parts[7]\t$parts[8]\n";
        $sid{$id}=1;
        if( $rcnt{$id} ) {
            $score{$id} += $parts[0];
            $rcnt{$id} += $parts[1];
        } else {
            $score{$id} = $parts[0];
            $rcnt{$id} = $parts[1];
        }
      }
      close(FILE);
   }
}

my %rank=readRank($rankFile);

foreach my $id (keys %sid) {
   my ($tid,$gid)=split(/ /,$id);
   my $pval=$rank{$gid};
   if( !$pval ) {
      print "error [$id]\n";
   }
   my $avg=$score{$id}/$rcnt{$id};
   print "$avg\t$rcnt{$id}\t$tid\t$pval\n";
}


sub readRank {
   my($file)=@_;
   my %res;
   my $zfile= new IO::Uncompress::Gunzip $file or die "gunzip failed: $GunzipError\n";
   while(my $line = $zfile->getline()) {
      chomp($line);
      my @vals=split(/\t/,$line);
      my $gid=$vals[1];
      $res{$gid}=$line;
   }
   close(FILE);
   return %res;
}
