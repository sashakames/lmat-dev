#!/usr/bin/perl -w

use strict;

my (%ans,%all_hdr);
my %fnd;
my ($tot_corr,$err_cnt) = (0,0);
my %rank_cnt;
my $rl_file=$ARGV[0];
my %allInfo=readAns($ARGV[1]);
my $type=$ARGV[2];
my $cutoff=$ARGV[3];
my (%incr_once,%corr_once);
my %wrong_rank_cnt;
my %merge_calls;
my ($no_hit_cnt,$low_conf_cnt,$uninf_cnt)=(0,0,0);
open(FILE,$rl_file) || die "fail $rl_file\n";
my $tot_uncl=0;
my %once;
while(my $line = <FILE>) {
   chomp($line);
   next if( $line !~ /\|s__(.*?)\t/ );
   my $node=$1; 
   $node =~ s/_/ /g;

   my $rcnt_sum = (split(/\t/,$line))[1];
   my $rval=$node;
   my $rank="species";
   if( $node =~ /Rickettsia conorii str. Malish/ ) {
      $rval="Rickettsia conorii";
      $rank=$1;
   } 
   elsif( $node =~ /Leptospira biflexa serovar Patoc/ ) {
      $rval="Leptospira biflexa";
      $rank=$1;
   } elsif( $node =~ /no_rank,Candidatus Baumannia/ ) { 
      $rval="Baumannia";
      $rank="genus";
   }
   $once{"$rank $rval"}=1;
   ## annoying tax tree quirk, classified in tree as 
   ## annoying synonyms that must be counted
   if( $rval =~ /Bacillus subtilis subsp/ ) {
      $rval = "Bacillus subtilis";
   }
   if( $rval =~ /Lactobacillus plantarum subsp./ ) {
      $rval = "Lactobacillus plantarum";
   }
   if( $rval =~ /Acidovorax citrulli/ ) {
      $rval = "Acidovorax avenae subsp. avenae";
   }
   if( $rval =~ /Borrelia bavariensis/ ) {
      $rval = "Borrelia garinii";
   }
   if( $rval eq "Candidatus Baumannia cicadellinicola") {
      $rval="Baumannia cicadellinicola";
   }
   if( ($rank eq $type) ) {
      $merge_calls{"$rval"} += $rcnt_sum;
      # record that parental values were classified
      #last;
   } 
}

foreach my $rval (keys %merge_calls) {
   my $no = $merge_calls{$rval};
   if($no<=$cutoff) {
      $tot_uncl += $no;
      next;
   }
   my $rank="";
   my $fndMe=0;
   foreach my $cand (keys %allInfo) {
      #if( $rval =~ /$cand/i ) {
      if( lc($rval) eq lc($cand)) {
         $fndMe=1;
         $fnd{$cand}=1;
      }
   }
   if( $fndMe ) {
      #print "lincorr:\t$rval\n";
      $fnd{$rval}=1;
      $rank_cnt{"$rval"} = "$no"; 
      $tot_corr++;
   } else {
      #print "linincr:\t$rval\n";
      $wrong_rank_cnt{"$rval"} = "$no"; 
      $err_cnt++;
   }
}
close(FILE);
my $miss_cnt=0;
foreach my $key (keys %allInfo) {
   if( !$fnd{$key} ) {
      print "missed: $key\t$allInfo{$key}\n";
      $miss_cnt++;
   }
}
my $fnd_cnt=0;
foreach my $key (keys %rank_cnt) {
      print "found:\t$key\t$rank_cnt{$key}\n";
      $fnd_cnt++;
}
my $wrng_cnt=0;
foreach my $key (keys %wrong_rank_cnt) {
      print "wrong:\t$key\t$wrong_rank_cnt{$key}\n";
      $wrng_cnt++;
}
print "final: corr: $fnd_cnt wrong: $wrng_cnt missed: $miss_cnt uncls: $tot_uncl\n";

sub isConsistent {
   my ($tname,$ri,$all) = @_;
   my $rank="";
   #print "check: [$tname] versus [$ri]\n";
   my @vals=split(/\t/,$ri);
   for(my $it = 1; $it <= $#vals; $it++) {
      if( $vals[$it] =~ /$tname/i ) {
         $all->{$vals[$it]} += 1;
         $rank=(split(/,/,$vals[$it]))[0];
         last;
      }
   }
   return $rank;
}

sub readAns {
   my($file)=@_;
   my %res;
   open(FILE,$file) || die "fail [$file]\n";
   while(my $line =<FILE>) {
      chomp($line);
      my $val=(split(/\s+/,$line))[0];
      $val =~ s/^[s|g]__//g;
      $val =~ s/_/ /g;
      $res{$val}=$line;
   }
   close(FILE);
   return (%res);
}
