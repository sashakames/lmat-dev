#!/usr/bin/perl -w 

use strict;


my $file_name=$ARGV[0];
procAll($file_name); #,$thresh,$cov_thresh);


sub procAll {
   my ($file) = @_;
   open(FILE,$file) || die "failed $file\n";
   open(OFILE,">$file.summary") || die "failed to write to $file.summary";
   <FILE>; ## skip first line
   while(my $line = <FILE>) {
      chomp($line);
      my @vals=split(/,/,$line);
      my $gcnt=int($vals[4]);
      next if( $gcnt < 1);
      my $spec="species,$vals[0] $vals[1]";
      my $tid=-1;
      my $par="";
      ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
      $vals[2] =~ s/\s+$//g;
      if( $tid == -1 ) {
         $spec="species,$vals[0] $vals[1]. $vals[2]";
         ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
      }
      if( $tid == -1 ) {
         $spec =~ s/\.//g;
         ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
      }
     if( $tid == -1 ) {
         $spec="no rank,$vals[0] $vals[1] $vals[2]";
         ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
      }
     if( $tid == -1 ) {
         $spec="no rank,$vals[0] $vals[1]. $vals[2]";
         ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
      }
     if( $tid == -1 ) {
         $spec="no rank,$vals[0] $vals[1] $vals[2]";
         $spec =~ s/_/ /g; 
         ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
      }
      my @ranks = ("species","no rank");
      for my $tr (@ranks) {
         if( $tid == -1 ) {
            $spec="$tr,$vals[0] $vals[1] $vals[2]";
            $spec =~ s/unfinished//g;
            $spec =~ s/chromosome//g;
            ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
         }
         if( $tid == -1 ) {
            $spec="$tr,$vals[0] $vals[1]. $vals[2]";
            ($tid,$spec,$par)=getLin($tid,$spec,$gcnt);
         }
      }
      my $lstr="$spec;$gcnt$par";
      print "[$tid][$spec][$gcnt][$par]\n";
      print OFILE "$tid\t$lstr\n";
   }
   close(OFILE);
   close(FILE);
}
   
sub getLin {
      my ($tid,$spec,$gcnt) = @_;
      my $cmd="grep \"$spec\" ~/repo/metag_repo/dev/runtime_inputs/ncbi_taxonomy_rank.txt |";
      my $lcal="";
      my $par="";
      open(P,$cmd) || die "fail $cmd\n";
      while(my $chk=<P>) {
         chomp($chk);
         my @tv=split(/\t/,$chk);
         $lcal=$tv[$#tv];
         if( $lcal eq $spec) {
            if( $tv[0] =~ /ktaxid=(\d+),/ ) {
               $tid=$1;
               my $sidx=$#tv-1;
               for(my $it=$sidx; $it > 0; $it--) {
                  $par = "$par:$tv[$it];0";
               }
            }
            last;
         }
         
      }
      close(P);
      #print "[$lcal] [$par]\n";
      my $pstr="$lcal-$gcnt-$par\n";
      #print "debug -> [tid=$tid] [$spec] [$par]\n";
      return ($tid,$spec,$par)
}
