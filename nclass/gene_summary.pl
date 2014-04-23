#!/usr/bin/perl -w 

use strict;

my $sz = @ARGV;
if ($sz != 5) {
  print "usage: echo <read-label-output-filename> | ./gene_summary.pl <rank-filename>  <log-odds-threshold> <min valid kmers for read> <min read length>\n";
  exit(1);
}


print "Read geneinfo: $ARGV[0]\n";
my %geneInfo=readGeneInfo($ARGV[0]);
print "Done geneinfo: $ARGV[0]\n";
my $sig_thresh=$ARGV[1];
my $min_kmers=$ARGV[2];
my $min_seq_len=$ARGV[3];
my $min_tax_score=$ARGV[4];


while(my $file_name=<STDIN>) {
   chomp($file_name);
   print "proc: $file_name\n";
   procAll($file_name,\%geneInfo,$sig_thresh,$min_kmers,$min_seq_len,$min_tax_score); 
}


sub procAll {
   my ($file,$rankRef,$thresh,$min_kmers,$min_seq_len,$min_tax_score) = @_;
   my %rank = %$rankRef;
   my %save_gene;
   my %save_tid;
   my %once;
   my $valid_read=1;
   my (%cnt_call,%child,%parent);
   open(FILE,$file) || die "failed $file\n";
   while(my $line = <FILE>) {
      chomp($line);
      my @vals=split(/\t/,$line);
      my ($tid,$tscore)=split(/ /,$vals[2]);
      my $lidx=$#vals;
      my ($ktaxid,$loscore,$label_type) = split(/ /,$vals[$lidx]);
      if( !$label_type ) {
         print "huh: $line\n";
         print "$vals[4]\n";
         print "$vals[3]\n";
         exit(0);
      }
      if( $label_type eq "ReadTooShort" ) {
         $cnt_call{"ShortRead"}++;
         next;
      }
      if($label_type eq "NoDbHits" ) {
         $cnt_call{"NoMatch"}++;
         next;
      }
      if($label_type eq "LCA_ERROR" ) {
         $cnt_call{"LCA_ERROR"}++;
         next;
      }
      
      if( $min_kmers != -1 ) {
         my $valid_kmers = (split(/ /,$vals[3]))[2];
         if( !$valid_kmers) {
            print "unexpected error [$vals[3]] [$line]\n";
            exit(0);
         }
         if( $valid_kmers < $min_kmers) {
            $cnt_call{"ShortRead"}++;
            next; 
         }
      }
      if( $min_seq_len != -1 ) {
         my $read_len = length($vals[1]);
         if( $read_len < $min_seq_len) {
            $cnt_call{"ShortRead"}++;
            next; 
         }
      }
      if($loscore < $thresh ) {
         $cnt_call{"LowScore"}++;
         next;
      }
      if($ktaxid <= 0 ) {
         print "When taxid=[$ktaxid] is <= 0, a previous category should have been triggered\n";
         print "$line\n";
         exit(0);
      }
      my $lineage=$rank{$ktaxid};
      if( !$lineage) {
         if( !$once{$ktaxid} ) {	
            print "here? $ktaxid\n";
            $once{$ktaxid} = 1; 
         }
         next;
      }
      my @lv=split(/\t/,$lineage);
      my $curr="root";
      if( $#lv == 0 ) {
         $curr="root"; 
      } else {
         $curr=$lineage;
      }
      $cnt_call{$curr}++;
      $save_gene{$curr} = $ktaxid;
      $save_tid{$curr} = $tid;
      if( $tscore >= $min_tax_score ) {
         $tax_cnt_call{$curr}++;
         $tax_save_gene{$curr} = $ktaxid;
         $tax_save_tid{$curr} = $tid;
      }
   }
   close(FILE);
   open(OFILE,">$file.$sig_thresh.$min_seq_len.genesummary") || die "failed to write to $file.summary";
   foreach my $node (keys %cnt_call) {
      my $node_cnt=$cnt_call{$node};
      my $lstr="";
      if( !$save_tid{$node} ) {
         $lstr="NULL\t$node;$node_cnt";
      } else {
         $lstr="$save_tid{$node}\t$node;$node_cnt";
      }
      print OFILE "$lstr\n";
   }
   close(OFILE);
   open(OFILE,">$file.$sig_thresh.$min_seq_len.genesummary.$min_tax_score.taxid") || die "failed to write to $file.summary";
   foreach my $node (keys %tax_cnt_call) {
      my $node_cnt=$tax_cnt_call{$node};
      my $lstr="";
      if( !$tax_save_tid{$node} ) {
         $lstr="NULL\t$node;$node_cnt";
      } else {
         $lstr="$tax_save_tid{$node}\t$node;$node_cnt";
      }
      print OFILE "$lstr\n";
   }
   close(OFILE);
}
   
sub readGeneInfo {
   my($file)=@_;
   my %res;
   open(FILE,$file) || die "fail $file\n";
   while(my $line =<FILE>) {
      chomp($line);
      next if( $line =~ /^#/ );
      my @vals=split(/\t/,$line);
      my $id=$vals[1];
      $res{$id}=$line; #"$vals[2]\t$vals[8]\t$vals[9]";
   }
   close(FILE);
   return %res;
}
