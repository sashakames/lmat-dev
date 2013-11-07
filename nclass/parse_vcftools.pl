#!/usr/bin/perl 

# example:    parse_bcftools.pl $ref $alignout.bcf  variants.$i.xml


my $fasta= $ARGV[0];
my $bcf=$ARGV[1];
my $out=$ARGV[2];

my @ids=();

open VARIANTS,">$out";
my %strain=&read_fasta_input($fasta);

$flank=50;
open IN,"$bcf";
while (my $line=<IN>) {
    if ($line !~ /^#/) {
	@cols=split/\t/,$line;
	$id=$cols[0];
	$pos=$cols[1];
	$ref=$cols[3];
	$var=$cols[4];
	$lv=length($var)-1;
	$lr=length($ref)-1;
	 
	if ($lv > $lr) { #$var=~ /\+/) {
	    $ref=substr($ref,1);
	    $var=substr($var,1);

	    $ref .="-" x $lv;
	    $context = substr($strain{$id},$pos-$flank,$flank).$var.substr($strain{$id},$pos+$lr,$flank);
	} elsif ($lv < $lr) { #($var=~/-/) {
	    $ref=substr($ref,1);
	    $var=substr($var,1);

	    $var .="-" x $lr;

	    $context = substr($strain{$id},$pos-$flank,$flank).substr($strain{$id},$pos+$lr,$flank);
	} else {
	    $context = substr($strain{$id},$pos-1-$flank,$flank).$var.substr($strain{$id},$pos,$flank);
	}
	$offset=$flank+1;
	print VARIANTS "<sequenceDifference>\n";
	print VARIANTS "<variant ref_name=\"$id\" ref_start=\"$pos\" ref_seq=\"$ref\">$var</variant>\n";
	print VARIANTS "<context offset=\"$offset\">$context</context>\n";
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

