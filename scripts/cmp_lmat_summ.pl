#!/usr/bin/perl -w

use strict;

my $f1=$ARGV[0];
my $f2=$ARGV[1];
my %rank=readRank($ARGV[2]);

my %all_specp_calls;
my %d1=loadSumm($f1,\%all_specp_calls,\%rank);
my %d2=loadSumm($f2,\%all_specp_calls,\%rank);

my ($agree_cnt,$d1_only,$d2_only)=(0,0,0);
foreach my $call (keys %all_specp_calls) {
   if( $d1{$call} && $d2{$call} ) {
      $agree_cnt++;
      print "agree $call\n";
   } elsif( $d1{$call} ) {
      print "d1 only $call\n";
      $d1_only++;
   } elsif( $d2{$call} ) {
      print "d2 only $call\n";
      $d2_only++;
   }
}

print "Result: $agree_cnt $d1_only $d2_only\n";
 
sub loadSumm {
   my ($file,$saveRef,$rankRef) = @_;
   my %save_calls;
   open(FILE,$file) || die "fail [$file]\n";
   while(my $line = <FILE>) {
      chomp($line);
      next if( $line=~ /^Abun/);
      my @vals=split(/\t/,$line);
      my $tid=$vals[7];
      my $call=$vals[8];
      if( $vals[8] =~ /no rank/) {
         ## get rank value call for now
         if( $rankRef->{$tid} =~ /\tspecies,(.*)\t/ ) {
            $call="species,$1";
         } elsif( $rankRef->{$tid} =~ /\tspecies,(.*)$/ ) {
            $call="species,$1";
         } elsif($rankRef->{$tid} =~ /\tgenus,(.*)\t/ ) {
            $call="genus,$1";
         } elsif( $rankRef->{$tid} =~ /\tgenus,(.*)$/ ) {
            $call="genus,$1";
         } else {
            print "Warning, really no rank is higher than genus?\n";
         }
      }
      $save_calls{$call} = $line;
      if( !$saveRef->{$call} ) {
         $saveRef->{$call} = $line;
      }
   }
   close(FILE);
   return %save_calls;
}


sub readRank {
   my($file)=@_;
   my %res;
   open(FILE,$file) || die "fail [$file]\n";
   while(my $line =<FILE>) {
      chomp($line);
      $line =~ s/no rank/no_rank/g;
      my @vals=split(/\t/,$line);
      my ($ntid,$ktid)=("","");
      if( $vals[0]=~ /,taxid=(\d+),ktaxid=(\d+),/) {
         ($ntid,$ktid)=($1,$2);
      }
      $res{$ntid}=$line;
   }
   close(FILE);
   return %res;
}

