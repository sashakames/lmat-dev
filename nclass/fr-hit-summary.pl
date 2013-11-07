#!/usr/bin/perl -w 

use strict;

my %rank=readRank($ARGV[0]);
#my $thresh=$ARGV[1];
#my $cov_thresh=$ARGV[1];

while(my $file_name=<STDIN>) {
   chomp($file_name);
   print "proc: $file_name\n";
   procAll($file_name,\%rank); #,$thresh,$cov_thresh);
}


sub procAll {
   my ($file,$rankRef) = @_;
   my %rank = %$rankRef;
   my %save_taxid;
   my %once;
   my (%cnt_call,%child,%parent);
   open(FILE,$file) || die "failed $file\n";
   while(my $line = <FILE>) {
      chomp($line);   
      my @vals=split(/\t/,$line);
      my $gtaxid=-1;
      if( $vals[1] =~ /\((\d+):/ ) {
         $gtaxid=$1;
      } else {
         if( $line !~ /UNKNOWN/ ) {
            print "unexpected error: $line\n";
         }
         next;
      }
      my $lineage=$rank{$gtaxid};
      if( !$lineage) {
         if( !$once{$gtaxid} ) {	
            print "huh: $gtaxid [$line]\n";
            $once{$gtaxid} = 1; 
         }
         next;
      }
      my @lv=split(/\t/,$lineage);
      for(my $it=1; $it <= $#lv; $it++) {
         my $curr=$lv[$it];
         if( $it == $#lv ) {
            $cnt_call{$curr}++;
            $save_taxid{$curr} = $gtaxid;
         }
         if( $it + 1 <= $#lv ) {
            my $cval=$lv[$it+1];
            $child{$curr} .= ":$cval";
         }
         if( $it > 1 ) {
            my $pval=$lv[$it-1];
            if( !$parent{$curr} ) {
               $parent{$curr} = $pval;
            } else {
               if( ! ($parent{$curr} eq $pval ) ) {
                  my $op=$parent{$curr};
                  print "unexpected tree [$curr] [$op] [$pval] [$line]\n";
               }
            }
         }
      }
   }
   close(FILE);

   open(OFILE,">$file.summary") || die "failed to write to $file.summary";
   foreach my $node (keys %cnt_call) {
      if( !$child{$node} ) {
         my $leaf_node_cnt=$cnt_call{$node};
         my $tot_read_cnt=$leaf_node_cnt;
         my $lstr="$save_taxid{$node}\t$node;$leaf_node_cnt";
         my $p_node=$parent{$node};
         while( $p_node ) {
            my $p_node_cnt=0;
            $p_node_cnt=$cnt_call{$p_node} if($cnt_call{$p_node});
            $tot_read_cnt += $p_node_cnt;
            $lstr .= ":$p_node;$p_node_cnt";

            $p_node = $parent{$p_node};
         }
         print OFILE "$lstr\n";
      }
   }
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
      my $tid="";
      if( $vals[0]=~ /,taxid=(\d+),/) {
         $tid=$1;
      }
      $res{$tid}=$line;
   }
   close(FILE);
   return %res;
}
