#!/usr/bin/perl -w

use strict;

my %add;
my $tmpfile="tmp.txt";
my $errlog="errlog.txt";
open(OFILE,">$tmpfile") || die "fail $tmpfile\n";
open(ERRLOG,">$errlog") || die "fail $errlog\n";
while(my $line = <STDIN>) {
	chomp($line);
	my ($gid,$tid,$kid) = (-1,-1,-1);
	if( $line =~ /gi\|(\d+)/ ) {
		$gid=$1;
	}
	if( $line =~ /\[tax_node_id (\d+)\]/ ) {
		$tid=$1;
        }	
	if( $line =~ /kp\|(\d+)/ ) {
		$kid=$1;
        }	
	if( ($kid == -1 || !$kid) && ($tid == -1 || !$tid) ) {
		print ERRLOG "error: $line\n";
	} else {
		if( $tid == -1 ) {
			if( $kid == -1  ) {
				print ERRLOG "error $line\n";
				exit(0);	
			}
			print OFILE "$kid\n";
			if(!$add{$kid} ) { 
				$add{$kid} = "$gid\t$kid\t$line";
			} else {
				$add{$kid} .= "\n$gid\t$kid\t$line";
			}	
	        } else {	
			print "$tid\t$gid\t$kid\t$line\n";
		}
	}
}
close(OFILE);
close(ERRLOG);
my $cmd = "./bulkLookupTaxNodeBySeqDataID.pl $tmpfile |";
open(PIPE,$cmd) || die "fail $cmd\n";
while(my $line = <PIPE> ){
	chomp($line);	
	my ($kid,$tid)=split(/\t/,$line);
	my @lines=split(/\n/,$add{$kid});
	foreach my $entry (@lines) {
		if(!$entry || !$tid) {	
			print "error: [$entry] [$tid] [$add{$kid}] [$kid]\n";
		 	exit(0);	
				
		}
		print "$tid\t$entry\n";
        } 
}
