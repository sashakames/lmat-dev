#!/usr/bin/perl -w 

use strict;

my %rank=readRank($ARGV[0]);
my $file_name=$ARGV[1];
my $ofile=$ARGV[2];
procAll($file_name,\%rank,$ofile);


sub procAll {
   my ($file,$rankRef,$ofile) = @_;
   my %rank = %$rankRef;
   my %save_taxid;
   my %once;
   my (%cnt_call,%child,%parent);
   my $read_cnt=0;
   open(OFILE,">$ofile") || die "failed $ofile\n";
   open(FILE,$file) || die "failed $file\n";
   while(my $line = <FILE>) {
      next if( $line !~ /^read_label/);
      chomp($line);   
      my @vals=split(/ /,$line);
      my $lidx=$#vals;
      my $ktaxid=$vals[$lidx-1];
      my $prob=$vals[$lidx];
      my $cov=$vals[2];
      my $ntid=-1;
      my $conf=$prob*$cov;
      if( $ktaxid != -1 ) {
          $ntid=$rank{$ktaxid};
          if( !$ntid ) {
            print "no mapping for $ktaxid\n";
            $ntid=-1; 
          }
      }
      print OFILE "read_$read_cnt\t$ntid\t$conf\n";
      $read_cnt++;

   }
   close(FILE);
   close(OFILE);
}
   
sub readRank {
   my($file)=@_;
   my %res;
   open(FILE,$file) || die "fail $file\n";
   while(my $line =<FILE>) {
      chomp($line);
      $line =~ s/no rank/no_rank/g;
      my @vals=split(/\t/,$line);
      my ($ntid,$ktid)=("","");
      if( $vals[0]=~ /,taxid=(\d+),ktaxid=(\d+),/) {
         ($ntid,$ktid)=($1,$2);
      }
      $res{$ktid}=$ntid;
   }
   close(FILE);
   return %res;
}
