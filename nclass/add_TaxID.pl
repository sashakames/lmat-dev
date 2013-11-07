#!/usr/bin/perl 

# ~/script/fasta_addTaxonomy.pl fasta.table /usr/mic/bio/shea/KPATH_TAXONOMY
# ~/script/fasta_addTaxonomy.pl fasta
# ~/script/fasta_addTaxonomy.pl blastout.table

# optional 3rd argument. default is /usr/mic/bio/shea/KPATH_TAXONOMY

#  ~/script/getTaxonomyFor_taxID.pl KPATH_TAXID > ! KPATH_TAXONOMY
#
## Rename sequences:
## On instdb-a: cdbin; kpathID_taxIDtable_simple.pl > ! KPATH_TAXID
## scp KPATH_TAXID to /usr/mic/bio/shea on yana; 
# cd /usr/mic/bio/shea; ~/script/getTaxonomyFor_taxID.pl KPATH_TAXID > ! KPATH_TAXONOMY ####
# ~/script/fasta.table_addTaxonomy.pl fasta.table /usr/mic/bio/shea/KPATH_TAXONOMY
## produces fasta.table.taxonomy
## Check if sequences have good names already. If they don't, and they start with kp|#seqID| then replace with taxonomy names
#set z=`grep -c '>kp|' fasta.table`
#echo "num seqs with kp names instead of genus,species,or strains: $z"
#if ( $z > 0) then
#    mv fasta.table fasta.table-original
#    mv fasta.table.taxonomy fasta.table
#endif

use Bio::LITE::Taxonomy::NCBI;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;

use Bio::DB::Taxonomy;

my $gi2taxid="/usr/mic/bio/shea/Taxonomy/gi_taxid_nucl.dmp";
my $namesfile="/usr/mic/bio/shea/Taxonomy/names.dmp";
my $nodesfile="/usr/mic/bio/shea/Taxonomy/nodes.dmp";


my $gi2taxid_binary="/usr/mic/bio/shea/Taxonomy/gi_taxid_nucl.bin";

#don't need, faster new way with Bio::LITE
#my $dbh = Bio::DB::Taxonomy->new(-source   => 'entrez'); 


# you need this for gi2taxid using Bio::LITE::Taxonomy::NCBI::Gi2taxid 
if (!-s $gi2taxid_binary) {
    new_dict (in => "$gi2taxid",
	      out => "$gi2taxid_binary");  
}

my $db = Bio::LITE::Taxonomy::NCBI->new (
					 db=>"NCBI",
					 names=> "$namesfile",
					 nodes=>"$nodesfile",
					 dict=>"$gi2taxid_binary"
					 );


my $in=$ARGV[0];  # fasta or list (eg blast outfmt 6 or 7, or fasta.table) with >kp|seqDataID|.... or gi|#
my $out=$in.".taxonomy";

my $kpath2taxonomyList;

if ($ARGV[1]) {
    $kpath2taxonomyList= $ARGV[1]; # /usr/mic/bio/shea/KPATH_TAXONOMY
} else {
    $kpath2taxonomyList="/usr/mic/bio/shea/KPATH_TAXONOMY";
}
my %kpSeqDataIDhash=();
my %kpSeqIDhash=();


open IN,"$kpath2taxonomyList";
while (my $line=<IN>) {
    if ($line =~ /seqDataID=(\d+)\t/) {
	$kpSeqDataIDhash{$1}=$line;
    }
    if ($line =~ /seqID=(\d+)\t/) {
	$kpSeqIDhash{$1}=$line;
    }
}

close IN;


open IN,"$in";
open OUT,">$out";

my $kpathSeqID;
my $kpathSeqDataID;
my $gi;
my $tax_node_id;
my $taxid="";

while (my $line=<IN>) {
    chomp $line;
    
    
    $kpathSeqID=$kpathSeqDataID=$tax_node_id="";
    if ($line =~ /kpath_id\|(\d+)/) {
	$kpathSeqDataID=$1;
    }
    
    if ($line =~ /kp\|(\d+)\|(\d+)/) {
	$kpathSeqID=$2;
	$kpathSeqDataID=$1;
    } 
    if ($line =~ /\[seq_data_id\s(\d+)\]/) {
	$kpathSeqDataID=$1;
    }
    if ($line =~ /\[sequence_id\s(\d+)\]/) {
	$kpathSeqID=$1;
    }
    if ($line =~ /gi\|(\d+)/ || $line =~ /gi(\d+)/ ) {
	$gi=$1;
    } else {
	$gi="";
    }
    
    if ($line =~ /\[tax_node_id\s(\d+)\]/) {
	$tax_node_id=$1;
    } 
    
    if (($kpathSeqID ne "" || $kpathSeqDataID ne "" || $tax_node_id  ne "" ) && $gi eq "")  {
	$count++;
    }
    
    #print "\$kpathSeqDataID=$kpathSeqDataID\n";
    if (defined $kpSeqDataIDhash{$kpathSeqDataID}) {
	$taxinfo=$kpSeqDataIDhash{$kpathSeqDataID};
	if ($taxinfo =~ /taxID=(.*?)\t/) {  
	    $taxid=$1;
	}
    } elsif (defined $kpSeqIDhash{$kpathSeqID}) {
	$taxinfo=$kpSeqIDhash{$kpathSeqID};
	if ($taxinfo =~ /taxID=(.*?)\t/) {  
	    $taxid=$1;
	}
    } elsif  ($gi ne "" || $tax_node_id ne "")   {
	if ($tax_node_id eq "" ) {
	    $taxid = $db->get_taxid($gi);
	} else {
	    $taxid = $tax_node_id;
	}
    }
    if ($line =~ /^(\S+)\s(.*)/) {
	my $part1=$1;
	my $part2=$2;
	print OUT "$taxid\t$gi\t$part1\t$part2\n";
    } 

}

