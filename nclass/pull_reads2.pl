#!/usr/bin/perl -w
use MyUtils;
use FileHandle;

my $ofilebase=$ARGV[0];
my $id_file=$ARGV[1];
my $src=$ARGV[2];
my $thresh=$ARGV[3];
my $min_kmer=$ARGV[4];
my $min_seq_len=$ARGV[5];
my $select_type="exclusive"; ## make this the default

if( !($select_type eq "exclusive") && !($select_type eq "inclusive") ){ 
   print "select type must be specified as exclusive or inclusive\n";
   exit(0);
}
my %ofhall;
my %valid;
open(FILE,$id_file) || die "fail $id_file\n";
while(my $line=<FILE>) {
   chomp($line);
   if ( $select_type eq "inclusive" ) {
      $check{$line}=1;
   } else {
      my @vals=split(/ /,$line);
      my $main_id=$vals[0];
      $valid{$main_id}=$main_id;
      for(my $it=1; $it <= $#vals; $it++) {
         $valid{$vals[$it]} = $main_id;
      }
      my $ofile="$ofilebase.$main_id"; 
      $ofhall{$main_id}= FileHandle::new();
      $ofhall{$main_id}->open(">$ofile") || die "failed to create $ofile\n";
   } 
}
close(FILE);
my $fh;
if( $select_type eq "inclusive" ) {
   $fh = FileHandle::new();
   $fh->open(">$ofile") || die "fail $ofile\n";
}
my $cnt=0;
while(my $line =<STDIN>) {
   chomp($line);
   my @vals=split(/\t/,$line);
   my $rlen=length($vals[1]);
   my ($tid,$score,$type)=split(/ /,$vals[4]);
   my ($ig1,$ig2,$valid_kmers) = split(/ /,$vals[2]);
   if( $valid{$tid} && $score >= $thresh && $valid_kmers >= $min_kmer && ($rlen >= $min_seq_len && !($vals[1] eq "X") ) )  {
         my $read=$vals[1];
         my $ofh = $ofhall{$valid{$tid}};
         my $sl=length($read);
         my $hdr="$vals[0];tid=$tid;score=$score;mtype=$type;valid_kmers=$valid_kmers;uid=$cnt;src=$src";
         $cnt+=1;
         writeFasta(\$read,$hdr,$ofh);
   }
}

if( $select_type eq "exclusive" ) {
   foreach my $key (keys %check) {
      $check{$key}->close();
   }
} else {
   $fh->close();
}
