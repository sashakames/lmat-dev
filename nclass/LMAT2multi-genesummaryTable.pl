#!/usr/bin/perl

# LMAT2multi-genesummaryTable.pl  -in flst  -out outfile  -min_reads 1 -min_frac 0.00 -min_score 0.7 
# Example:
# LMAT2multi-genesummaryTable.pl /p/lscratchf/shea/Gene2Megan/genesummary.list..male_Tongue_dorsum -out male_Tongue_dorsum.megan   -min_reads 1 -min_frac 0.00 -min_score 0.5 


use Getopt::Long; 

my $min_reads=1;  # minimum number of raw reads to count in output
my $min_frac=0.000; # minimum fraction of total reads to count in output
my $file_list;  # file with two columns, listing the full path to each LMAT output file( .fastsummary or .species or .genus or .plasmid files), one path per line.
my $out;  # table summarizing output of files in -in file. taxID's are rows, 
my $avg_score_min=0.5;

GetOptions ( "min_reads:i" =>\$min_reads,
             "min_frac:f" =>\$min_frac,
             "min_score:f" =>\$avg_score_min,
             "in:s" =>\$file_list,
             "out:s" => \$out
    );
        
if($file_list !~ /\w/ || $out  !~ /\w/) {
    die "\nUsage: $0 -in fileList -out outfile \n\t-min_score <minimum average LMAT score per read for the tax id to be tabulated, default 0.5>\n\t-min_frac <minimum fraction of total reads in sample for a tax id to be tabulated, default 0>\n\t-min_reads <minimum number reads in sample for a tax id to be tabulated, default 1>\n\nExample: LMAT2multi-genesummaryTable.pl  -in fileList  -out genes.megan -min_frac 0 -min_score 0.5 -min_reads 1\n\nfileList is a file listing the files to include in the output table, where each line of fileList has two tab-delimited columns: 1) sample name and 2) full path name to the LMAT .genesummary file for that sample\n
KEGG and SEED ID's are in the rows and samples are in the columns. This program requires the LMAT *.log file which gives the total number of reads for that sample to be in the same directory as the LMAT .genesummary file.

The first time the script is used, first please download the following data files to the $LMAT_DIR
wget ftp://gdo-bioinformatics.ucllnl.org/lmat/gi2kegg.map.gz
wget ftp://gdo-bioinformatics.ucllnl.org/lmat/gi2seed.map.gz
wget ftp://gdo-bioinformatics.ucllnl.org/pub/lmat/geneID2proteinGI.gz
\n";
}

my $dir=$ENV{LMAT_DIR};
my $proteinGI2kegg="gi2kegg.map";
my $proteinGI2seed="gi2seed.map";
my $geneID2proteinGI="geneID2proteinGI";
if( $dir ) {
   $proteinGI2kegg="$dir/gi2kegg.map";
   $proteinGI2seed="$dir/gi2seed.map";
   $geneID2proteinGI="$dir/geneID2proteinGI";
}

my %kegg=();
my %seed=();
open IN,"$proteinGI2kegg" or die "Cannot open $proteinGI2kegg: $!.\n";
my $kegg_cnt=0;
while (my $line=<IN>){
    chomp $line;
    my @x=split/\t/,$line;
    my $proteinGI=$x[0];
    $proteinGI =~ s/GI//;
    my $keggID=$x[1];
    $kegg{$proteinGI}=$keggID;
    $kegg_cnt++;
}
close IN;
my $seed=0;
open IN,"$proteinGI2seed" or die "Cannot open $proteinGI2seed: $!.\n";
while (my $line=<IN>){
    chomp $line;
    my @x=split/\t/,$line;
    my $proteinGI=$x[0];
    my $seedID=$x[1];
    $seed{$proteinGI}=$seedID;
    $seed_cnt++;
}
close IN;


open IN,"$geneID2proteinGI";
my %geneID2proteinGI=();
my $id_cnt=0;
while  ( my $line=<IN>) {
    chomp $line;
    my ($geneID,$proteinGI)=split/\s/,$line;
    $geneID2proteinGI{$geneID}=$proteinGI;
   $id_cnt++;
}
close IN;

print "kegg_cnt=$kegg_cnt ; seed_cnt=$seed_cnt ; gid2pgi=$id_cnt\n";

open IN,"$file_list" or die "Cannot open $file_list: $!\n";
my $count=0;
my @File=();
my @Name=();
while ( my $line=<IN>) {
    chomp $line;
    ($name,$file)=split/\t/,$line;
    $Name[$count]=$name;
    $File[$count]=$file;
    $count++;
    print "$name\t$file\n";
}

close IN or warn "Cannot close $file_list: $!\n";

my $total_reads=0;
my %data_by_org=();

my @num_reads_array=();

for (my $i=0 ; $i<@File ; $i++) {
    my $file=$File[$i];

    $file =~ s/\*//;
    print "file: $file\n";
    my $file_nopath=`basename $file`;
    chomp $file_nopath;
    print "file_nopath: $file_nopath\n";
    my $prefix;
    if ($file_nopath =~ /^(\S+)_output/){
	$prefix=$1; 
    }
    my $thisdir=`dirname $file`;
    chomp $thisdir;
    print "dir: $thisdir\n";
    my $cmd="grep \" reads in\" $thisdir/*_output.log "  ;
    #print "$cmd\n";
    my $z=`$cmd`;
    my $num_reads=0;
    if ($z=~ /(\d+)\sreads\sin/) {
	$num_reads=$1;
    }
    print "num reads: $num_reads\n";
    $num_reads_array[$i]=$num_reads;
    $total_reads +=$num_reads;
    open IN,"$file" or warn "Cannot open $file: $!\n";
    my $count_gene_classified_reads=0;
    my $count_reads_classified_by_seed=0;
    my $count_reads_classified_by_kegg=0;
    my ($avg_score,$wrc,$rc,$taxid,$org);
    while (my $line=<IN>) {
	chomp $line;

#*.genesummary - read count assigned to gene, LMAT NCBI taxonomic ID, NCBI taxonomic ID for source gene, NCBI Gene ID, NCBI Locus tag, description, gene type and protein accession.
	
	
#	($rc,$taxid,$NCBI_taxid,$geneID,$locus_tag,$description,$gene_type,$protein_accession) = split/\t/,$line;
	($avg_score,$rc,$taxid,$NCBI_taxid,$geneID,$locus_tag,$description,$gene_type,$protein_accession) = split/\t/,$line;
	my $proteinGI=$geneID2proteinGI{$geneID};
	my $keggID=$kegg{$proteinGI};
	my $seedID=$seed{$proteinGI};
#	print "\$proteinGI=$proteinGI\tkeggID=$keggID\tseedID=$seedID\n";
#	print "avg_score=$avg_score\t\$rc=$rc\ttaxid=$taxid\tncbi_taxid=$NCBI_taxid\tgeneID=$geneID\tlocus_tag=$locus_tag\tdescription=$description\tgene_type=$gene_type\tprotein_accessio=$protein_accession\n";
	if ($rc/$num_reads >= $min_frac && $rc >= $min_reads  ) {
	    $count_gene_classified_reads += $rc;
	    $data_by_SEED{$seedID}{$i}=$rc;
	    $data_by_KEGG{$keggID}{$i}=$rc;
	    if ($seedID =~ /\d/) {
		$count_reads_classified_by_seed += $rc;
	    }
	    if ($keggID =~ /\d/) {
		$count_reads_classified_by_kegg += $rc;
	    }
	}
    }
    print "$dir\t$i\tLMAT gene classified=$count_gene_classified_reads\tSEED classified=$count_reads_classified_by_seed\tKEGG classified reads=$count_reads_classified_by_kegg\n";
    close IN;
} # for (my $i=0 ; $i<@File ; $i++) {
 
open OUT,">$out";
if ($out =~ /megan/i) {
    print OUT "\@Creator\tLMAT\n";
    print OUT "\@CreationDate\t",`date`;
    print OUT "\@ContentType\tSummary4\n";
    print OUT "\@Names\t",join("\t",@Name),"\n";
    #print OUT "\@Uids\t@sra_list\n";
    print OUT "\@Sizes\t",join("\t",@num_reads_array),"\n";
    print OUT "\@TotalReads\t$total_reads\n";
    #print OUT "\@Collapse\n";
    print OUT "\@Algorithm\tSEED\tlmat\n";
    print OUT "\@Algorithm\tKEGG\tlmat\n";

    #print OUT "\@NodeStyle\t
} else {
     print OUT "DB\tID\t",join("\t",@Name),"\n";
}
delete $data_by_KEGG{""};
delete $data_by_SEED{""};

foreach my $keggID (sort {$a <=> $b} keys %data_by_KEGG) {
    print OUT "KEGG\t$keggID";
   
    for (my $i=0 ; $i<@File ; $i++) {
	if (!defined $data_by_KEGG{$keggID}{$i}  ) {
	    $data_by_KEGG{$keggID}{$i}=0;
	}
	print OUT "\t$data_by_KEGG{$keggID}{$i}";
    }
    print OUT "\n";
}
foreach my $seedID (sort {$a <=> $b} keys %data_by_SEED) {
    print OUT "SEED\t$seedID";
   
    for (my $i=0 ; $i<@File ; $i++) {
	if (!defined $data_by_SEED{$seedID}{$i}  ) {
	    $data_by_SEED{$seedID}{$i}=0;
	}
	print OUT "\t$data_by_SEED{$seedID}{$i}";
    }
    print OUT "\n";
}
close OUT;


