#!/usr/bin/perl
#
$fnd=$ARGV[0];
open(FILE,$fnd);
while(<FILE>) {
   chomp;
   my ($id,$pcnt)=split;
   $line=$_;
   $line =~ s/^\d+ //g;
   $save{$id} = $line;

}
close(FILE);
$all=$ARGV[1];
open(FILE,$all);
while(<FILE>) {
   chomp;
   my ($id,$cnt)=split;
   push(@allcnt,"$id $cnt");
}
my $min_sample_size = $ARGV[2];
#for($it=$#allcnt; $it >= 0; $it--) {
my $it = $#allcnt;
while($it >= 0) {
   my $sample_val="";
   my $sample_size = 0;
   for($j = $it; $j >= 0; $j--) {
      my ($id,$cnt)=split(/ /,$allcnt[$j]);
      if( $save{$id} ) { 
         my @vals=split(/ /,$save{$id});
         $sample_size += $#vals+1;
         if( !$sample_val) {
             $sample_val = "$save{$id}";
         } else {
            $sample_val .= " $save{$id}";
         }
      }
      #print "huh: $sample_size $min_sample_size\n";
      if($sample_size >= $min_sample_size || $j == 0) {
         $stop_pos = $j;
         last;
      }   
   }
   if( !$sample_val ) {
      print "what happened: $it,$stop_pos,$id,$cnt,$sample_size\n";
      exit(0);
   }  
   my $diff=0;
   my ($mid1,$mcnt1)=split(/ /,$allcnt[$it]);
   my ($mid2,$mcnt2)=split(/ /,$allcnt[$stop_pos]);
   my $iter_diff=$it-$stop_pos;
   my $kcnt_diff=$mcnt2-$mcnt1;
   my @vals = split(/ /,$sample_val);
   my @svals=sort { $a cmp $b } @vals;
   
   my $sum = 0;
   foreach my $val (@vals) {
      $sum += $val;
   }
   my $rsum = 0;
   my $mean = $sum / ($#vals + 1 );
   foreach my $val (@vals) {
      my $diff = $val - $mean;
      $rsum += ($diff ** 2);
   }
   $stdev =  ($rsum / $#vals) ** 0.5;
   my ($min,$max)  = ($svals[0],$svals[$#svals]);
   for($j = $it; $j >= $stop_pos; $j--) {
      my ($mid1,$mcnt1)=split(/ /,$allcnt[$j]);
      print "$mid1 $min $max $mean $stdev $iter_diff $kcnt_diff $sample_size\n";
   }
   

   $it = $stop_pos -1;
}
close(FILE); 

