#!/usr/local/bin/perl -w

################

use strict;
use warnings;
use Carp;

use CGI ":standard";

use CBNP::Constants ":standard", ":kpath";
use CBNP::Boilerplate ":standard";
use CBNP::Kpath::SequenceCache;

use CBNP::MetaLite::Oracle;

{
	#my $seqDataID = $ARGV[0];
    #my $seqDataID = cgiParamRequired("seq_data_id");

    my $metaLite = CBNP::MetaLite::Oracle->new();
    my $kpathHandle = $metaLite->getROHandle("kpath",
      {
        LongReadLen => $KPATH_TEXT_SIZE,
      });


  my $listFileName = shift;
  if(!defined $listFileName)
  {
    confess("list file required as first parameter");
  }

  my $idList_r = loadListFile($listFileName);
  
  my $tempSQL = "SELECT sequence_id FROM seq_data WHERE seq_data_id = ?";
  my $sequenceSTH = $kpathHandle->prepare($tempSQL);

  $tempSQL = "SELECT tax_node_id FROM sequence WHERE sequence_id = ?";
  my $nodeSTH = $kpathHandle->prepare($tempSQL);

  foreach my $currSeqDataID(@{$idList_r})
  {
    print "$currSeqDataID\t";

    eval
    {
      $sequenceSTH->execute($currSeqDataID);
      my ($currSequenceID) = $sequenceSTH->fetchall_arrayref_exists()->[0]->[0];

	  $nodeSTH->execute($currSequenceID);
	  my ($currTaxNodeID) = $nodeSTH->fetchall_arrayref_exists()->[0]->[0];
	
	  print "$currTaxNodeID\n";
    };
    if ($@)
    {
	    print "NOT_FOUND\n";
    }
  }
}

sub loadListFile
{
  my ($fileName) = @_;

  open(INLIST, "< $fileName") || 
    confess("failed to open list file \"$fileName\" for reading: $!");
  my @wholeListFile = <INLIST>;
  close(INLIST);

  my @listFile = ();
  foreach my $line(@wholeListFile)
  {
    chomp($line);
    if($line =~ /^\s*\#/ ||
       $line =~ /^\s*$/)
    {
      next;
    }

    push(@listFile, $line);
  }
  
  return \@listFile; 
}

