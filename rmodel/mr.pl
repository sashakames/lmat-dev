#!/usr/bin/perl -w
use strict;
my $rl=$ARGV[0];
my $num_reads=$ARGV[1];
my $verbose=0;
if($verbose) {
   open(DEBUG,">debug.txt");
}

my $gc_partitions=4;
my $partition_size=$num_reads/$gc_partitions;
my @range=();
for(my $it = 0; $it < $gc_partitions; $it++) {
   my $beg=$it;
   $beg /= $gc_partitions;
   $beg *= 100.0;
   my $end=$it+1;
   $end /= $gc_partitions;
   $end *= 100.0;
   push(@range,"$beg $end");
}
my @gcbase=("c","g");
my @atbase=("a","t");
for(my $ri = 0; $ri < $num_reads; $ri++) {
   my $parti = int($ri/$partition_size);
   my ($gc_beg,$gc_end)=split(/ /,$range[$parti]);
   my $gc_range=($gc_end-$gc_beg);
   my $gc_draw=$gc_beg+int(rand($gc_range)); 
   my $gc_pcnt=$gc_draw/100.0;
   my $gc_cnt = int($gc_pcnt*$rl);
   my @bases=();
   print DEBUG "gc_range = $gc_range, $gc_pcnt, $gc_cnt $parti $gc_beg $gc_end\n" if($verbose);
   for(my $i = 0; $i < $gc_cnt; $i++) {
      my $idx = int(rand(2));
      push(@bases,$gcbase[$idx]);
   }
   my $at_cnt = $rl-$gc_cnt;
   for(my $i = 0; $i < $at_cnt; $i++) {
      my $idx = int(rand(2));
      push(@bases,$atbase[$idx]);
   }
   shuffle(\@bases);
   my $read = "@bases";
   $read =~ s/ //g;
   #my $slen = length($read);
   #my %kmer_cnt;
   #for(my $it=$ksize; $it < $slen; $it++) {
   print ">read_$ri\n";
   print "$read\n";
}
close(DEBUG) if($verbose);

sub shuffle {
  my $array = shift;
  my $i;
  for ($i = @$array; --$i; ) {
    my $j = int rand ($i+1);
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }
}

