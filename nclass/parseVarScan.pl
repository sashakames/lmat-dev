#!/usr/bin/perl 

# example:    parseVarScan.pl $ref $alignout.2cns  variants.$i.xml


my $fasta= $ARGV[0];
my $cns=$ARGV[1];
my $out=$ARGV[2];

my @ids=();

open VARIANTS,">$out";
my %strain=&read_fasta_input($fasta);

open IN,"$cns";
while (my $line=<IN>) {
    if ($line !~ /Chrom\s+Position\s+Ref/) {
	@cols=split/\t/,$line;
	$id=$cols[0];
	$pos=$cols[1];
	$ref=$cols[2];
	$var=$cols[3];
	if ($var=~ /\+/) {
	  
	    $var=~ s/\+//;
	    $ref="-" x length($var);
	    $context = substr($strain{$id},$pos-15,15).$var.substr($strain{$id},$pos,15);
	} elsif ($var=~/-/) {
	    $var=~ s/-//;
	    my $l=length($var);
	    $ref=$var;
	    $var="-" x $l;

	    $context = substr($strain{$id},$pos-15,15).substr($strain{$id},$pos+$l,15);
	} else {
	    $context = substr($strain{$id},$pos-1-15,15).$var.substr($strain{$id},$pos,15);
	}
	$offset=14;
	print VARIANTS "<sequenceDifference>\n";
	print VARIANTS "<variant ref_name=\"$id\" ref_start=\"$pos\" ref_seq=\"$ref\">$var</variant>\n";
	print VARIANTS "<context offset=$offset>$context</context>\n";
	print VARIANTS "</sequenceDifference>\n";

	#print "$id  $pos  $ref  $var  $context $pos_rel_context\n";
    }
}
close IN;
close VARIANTS;
sub read_fasta_input {
    my $infile=shift;
    open INDATA,"$infile" or die "Can't open $infile: $!\n";
        
    #get all sequences and make a hash called strain with id and sequence info
    my $sequence = "";
    my $id;
    my $num_genomes=0;
    my %strain=();
    while (my $line = <INDATA>) {
        if ($line =~ /^>(.*)/) {
            if ($sequence ne "") {
                $ids[$num_genomes]=$id;
                $strain{$id}=uc($sequence); 
                $num_genomes++;
            }
            $id = $1;
            $id =~ s/\r//g; 
            $id =~ s/\t/ /g; 
            $id =~ s/\s+/ /g;
            chomp($id);
            $sequence = "";
        } else {
            chomp($line);
            $line =~ s/[\n\t\s\r]//g; 
            $sequence .= $line;
        }
        
    }   
    
    # Add last one
    if (defined $id) {
        $ids[$num_genomes]=$id;
        $strain{$id}=uc($sequence);
        $num_genomes++;
    }
    
    close INDATA or warn  "Can't close $infile: $!\n";
    return %strain;

} # end sub read_fasta_input 

