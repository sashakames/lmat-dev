#!/usr/bin/perl -w

use strict;

while(my $line = <STDIN>) {
    if($line =~ /<sequencer method=\"(.*?)\">(.*?)<\/sequencer>/) {
      print "$1 $2\n";
    }
}
