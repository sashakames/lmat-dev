#!/usr/bin/perl -w

use strict;
use FileHandle;

print "entering";

my $tax_call_file=$ARGV[0];
my $gene_call_file=$ARGV[1];
my $snp_call_file=$ARGV[2];
my $labeled_read_file=$ARGV[3];
my $ofile=$ARGV[4];
my $remap_id=$ARGV[5];
my $xml_header=$ARGV[6];
my $custom=$ARGV[7];
my $min_gene_reads=$ARGV[8];
my $read_los=$ARGV[9];
my %remap=loadRemap($remap_id);
my %reads=loadReads($labeled_read_file,\%remap,$read_los);
my $fh = FileHandle::new();

print "output file: $ofile";

$fh->open(">$ofile") || die "fail to write to [$ofile]\n";

if($xml_header eq "none") {
   print $fh "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
   print $fh "<organismsReport>\n";
   print $fh "   <dataset>\n";
   print $fh "      <datasetName>UNKNOWN</datasetName>\n";
   print $fh "   </dataset>\n";
} else {
   open(DFILE,$xml_header) || die "failed to open $xml_header\n";
   while(my $line = <DFILE>) {
      last if( $line =~ /<\/organismsReport>/ );
      print $fh $line;
   }
   close(DFILE);
}
my %mapids;
open(FILE,$custom) || die "failed to open custom mappings: $custom\n";
while(my $line = <FILE>) {
   chomp($line);
   my @vals=split(/\t/,$line);
   $mapids{$vals[0]}=$vals[1];
}
close(FILE);
print $fh "   <organisms>\n"; 
my ($tax_calls_ref,$read_sum,$human_read_sum)=loadTaxCalls($tax_call_file);
my %snp_calls=loadSnpCalls($snp_call_file);
my %gene_calls=loadGeneCalls($gene_call_file,\%remap,$min_gene_reads,$read_sum);
my @tax_calls=@$tax_calls_ref;
foreach my $call (@tax_calls) {
   my $taxid=(split(/\t/,$call))[1];
   next if($taxid <= 2 || $taxid == 131567 || $taxid == 2759);
   my $read_cnt=(split(/\t/,$call))[0];
   if( $call =~ /Homo sapiens/) {
      $read_cnt = $human_read_sum;
   }
   if( $taxid != 9606 && $call =~ /Homo sapiens/ ) {
      next;
   } else {
      my $pid=$taxid;
      if( $mapids{$taxid} ) {
         $pid=$mapids{$taxid};
      } 
      #print "gc: $taxid $gene_calls{$taxid}\n" if($gene_calls{$taxid});
      reportCal($call,$reads{$taxid},$gene_calls{$taxid},$snp_calls{$taxid},$read_sum,$read_cnt,$fh,$pid);
   }
}

print $fh "   </organisms>\n"; 
print $fh "</organismsReport>\n";


sub getRankInfo {
   my ($tax_call) = @_;
   my @vals=split(/\t/,$tax_call);
   my $sl=$#vals;
   my ($strain_val,$species_val,$genus_val) = ("","","");
   $vals[$sl] =~ s/plasmid //g;


   if( $vals[2] eq "strain" ) {
       ## this case the genus and species names are appended after the lineage
      my @t = split(/\s+/,$vals[$sl-2]);
      $strain_val = $t[$#t];
      ## need to deal with species group here
      $species_val = $vals[$sl-1];
      $genus_val = $vals[$sl];
   } elsif( $vals[2] eq "species" ) {
      my @t = split(/\s+/,$vals[$sl]);
      $species_val = $t[$#t];
      $genus_val = $vals[$sl-1];
   } elsif( $vals[2] eq "subspecies" ) {
      my @t = split(/\s+/,$vals[$sl-1]);
      $species_val = $t[$#t];
      $genus_val = $vals[$sl-1];
   } elsif( $vals[2] eq "genus" ) {
      $genus_val = $vals[$sl];
   }
   return ($genus_val,$species_val,$strain_val);
}

sub reportCal {
   my ($tax_call,$read_calls,$gene_calls,$snp_calls,$read_sum,$read_cnt,$fh,$prn_tid) = @_;
   my @vals=split(/\t/,$tax_call);
   my $read_pcnt = 0;
   $read_pcnt = $read_cnt/$read_sum if($read_sum>0);
   $read_pcnt *= 100.0;
   my $read_pcnt_str=sprintf("%.3f",$read_pcnt);
   my $taxid=$vals[1];


   my $lin_str="";

   my $idx_offset = 0;

   ## strains are tricky so we put 

   my $idx_max = $#vals;
   
   
   for(my $it=4; $it <= ($idx_max); $it++) {
      $lin_str .="$vals[$it]; " if($it < $idx_max);
      $lin_str .="$vals[$it]." if($it == $idx_max);
   }


   my ($genes,$variants,$reads)=("","",""); 

   print $fh    "      <organism>\n";
   print $fh    "         <relativeAmount count=\"$read_cnt\">$read_pcnt_str</relativeAmount>\n";
   print $fh    "         <taxonomy taxon_id=\"$prn_tid\">$lin_str</taxonomy>\n";

   my $call=$vals[$idx_max];

   my  ($strain, $species, $genus) = ("", "", "");

   print "$vals[3]\n";
   
   if ($tax_call =~ /Homo sapien/) {
       #print $fh "         <commonName>human</commonName>\n";       
       #print $fh "         <hostGenus>Homo</hostGenus>\n";
       #print $fh "         <hostSpecies>sapiens</hostSpecies>\n";
   }
   elsif( $tax_call =~ /\t(genus|species|strain)\t/ ) {
      ($strain, $species, $genus) = split(',',$vals[3]) ;
      print $fh "         <organismName>$call</organismName>\n";
      print $fh "         <genus>$genus</genus>\n" if( $genus ne '-' );
      print $fh "         <species>$species</species>\n" if( $species ne '-' );
      print $fh "         <strain>$strain</strain>\n" if( $strain ne '-' );
   } else {
      print $fh "         <nearestNeighbor>$call</nearestNeighbor>\n";
   }
   if( $gene_calls ) {
      print $fh "         <genes>\n";
      my @gvals=split(/\n/,$gene_calls);
      foreach my $gene (@gvals) {
         next if($gene =~ /NULL/);
         #my ($cnt_tid,$gid,$sname,$product,$ignore)=split(/\t/,$gene);
         my ($cnt_tid,$read_tid,$gene_tid,$gid,$lctag,$product,$type,$acc)=split(/\t/,$gene);
      #print $fh "            <gene protein_id=\"undef\" locus_tag=\"undef\" product=\"$product\" ref_name=\"$gid\">$sname</gene>\n";
      #print $fh "            <gene product=\"$product\" ref_name=\"$gid\">$sname</gene>\n";
      ## this should be for non-protein coding genes
      if( !$acc ) {
         $acc=$gid;
      }
      print $fh "            <gene product=\"$product\" ref_name=\"$acc\" locus_tag=\"$lctag\">$lctag</gene>\n";
      }
      print $fh "         </genes>\n";
   }
   #if( $snp_calls ) {
      #print $fh "         <variants>\n";
      #open(VFILE,$snp_calls) || die "failed to open [$snp_calls]\n";
      #while(my $vl=<VFILE>) {
         #print $fh "         $vl";
      #}
      #close(VFILE);
      #print $fh "         </variants>\n";
   #} 
   if( $read_calls && $tax_call !~ /Homo sapien/ ) { # && ($genus || $species || $strain)) {
      print $fh "         <reads>\n";
      my @rvals=split(/\t/,$read_calls);
      foreach my $rval (@rvals) {
      print $fh "            <sequence>$rval</sequence>\n";
      }
      print $fh "         </reads>\n";
   }
   print $fh "      </organism>\n";
      
}
 
sub loadTaxCalls {
   my ($file)= @_;
   open(FILE,$file) || die "loadTaxCall fail $file\n";
   my $sum = 0;
   my $hsum = 0;
   my @save;
   while(my $line = <FILE>) {
      chomp($line);
      next if( $line =~ /Homo sapiens/ );
      my @vals=split(/\t/,$line);
      if( $line =~ /Homo sapiens/ ) {
         $hsum += $vals[0];
      }
      $sum += $vals[0];
      push(@save,$line);
   }
   close(FILE);
   return (\@save,$sum,$hsum);
}



sub loadRemap {
   my ($file)= @_;
   open(FILE,$file) || die "missing: loadRemap $file\n";
   my %save;
   while(my $line = <FILE>) {
      chomp($line);
      my @vals=split(/\s+/,$line);
      for(my $it=1; $it <= $#vals; $it++) {
         $save{$vals[$it]} = $vals[0];
      }
   }
   return %save;
} 

sub loadReads {
   my ($file_lst,$remap_ref,$min_los) = @_;
   my %save;
   my %remap=%$remap_ref;
   open(FILE_LST,$file_lst) || die "loadReads failed $file_lst\n";
   my %nodups;
   my $miss=0;
   while(my $file= <FILE_LST>) { 
      chomp($file);
      open(FILE,$file) || die "loadReads1 failed $file\n";
      while(my $line = <FILE>) {
         chomp($line);
         my @vals=split(/\t/,$line);
         my $hdr=$vals[0];
         next if( !$hdr );
         my ($taxid,$los,$type)=split(/ /,$vals[$#vals]);
         next if( $los < $min_los || $taxid == 9606 ); ## don't bother with human
         if($taxid == 754086 ) {
            print "yes? $line\n";
         } 
         if( $nodups{$hdr} ) {
            my $pr=$remap{$taxid};
            $pr="none" if(!$remap{$taxid});
            print "warning found a duplicate, [$hdr] [$taxid] [$pr] [$nodups{$hdr}]\n";
            print "debug: $line\n";
            next;
         }
         if( $remap{$taxid} ) {
            $taxid=$remap{$taxid};
         }
         $nodups{$hdr}=$taxid;
         if( $save{$taxid} ) {
            $save{$taxid} .= "\t$hdr";
         } else {
            $save{$taxid} = "$hdr";
         }
      }
      close(FILE);
   }
   close(FILE_LST);
   return %save;
}

sub loadGeneCalls {
   my ($file,$remap_ref,$min_reads,$read_tot)= @_;
   my %remap=%$remap_ref;
   my %save;
   open(FILE,$file) || die "fail $file\n";
   while(my $line = <FILE>) {
      chomp($line);
      my @vals=split(/\s+/,$line);
      my $num_reads=$vals[0];
      my $pcnts=$num_reads / $read_tot;
      next if( $pcnts < $min_reads);
      my $taxid=$vals[1];
      #print "genedebug: $pcnts $min_reads $num_reads $read_tot\n";
      $taxid = $remap{$taxid} if( $remap{$taxid} );
      if( !$save{$taxid} ) {
         $save{$taxid} = $line;
      } else { 
         $save{$taxid} .= "\n$line";
      }
   }
   close(FILE);
   return %save;
}

sub loadSnpCalls {
   my ($file) = @_;
   my %file_lst;
   if(open(FILE,$file)) { 
      while(my $line = <FILE>) {
         chomp($line);
         if( $line =~ /variants\.(\d+)\./ ) {
            my $taxid=$1; 
            $file_lst{$taxid} = $line;
         } else {
            print "parse error $line failed to get taxid\n";
         }
      } 
   } else {
      print "warning no snp file: $file\n";
   } 
   return %file_lst;
}
