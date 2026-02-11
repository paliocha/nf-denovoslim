#!/usr/bin/perl

use Getopt::Long;
use strict;

my $in;
my $out;
my $help;
my $helpMsg="
*********************************************
fastafix.pl
This script removes empty lines from fasta
file and writes the sequences on single rows.
Use: perl fastafix.pl -i str (-o str)
params:
-i string    Input fasta file, reguired.
-o string    Output fasta file, must not
             exist yet. Default output file: 
             fixed_inputfilename. Optional.
-h           Print help message and exit.
*********************************************
";
GetOptions ('i:s' => \$in,'o:s' => \$out, 'h' => \$help );
if($help==1) { print $helpMsg; exit; }
open(IN,$in) or die ("Give input fasta file as param -i");

#read the first line and check that it is a fasta header
$_=<IN>; 
if($_!~m/^>/) { die("Input file doesn't begin with fasta heading"); }

#output file open
if(length($out)>0) {
	if(-e $out) { die("Can't overwrite existing output file"); }
	open (OUT,">$out") or die("Can't open output file $out");
}#endif
#if output file parameter isn't given then use default output file value (fixed_$in)
else {
	if(-e "fixed_$in") { die("Can't overwrite existing default output file fixed_$in"); }
	open(OUT,">fixed_$in") or die("Can't open default output file fixed_$in");
}

#read other lines
print OUT $_; #print first header without nextrow 
while(<IN>) {
	if( m/^>/) {
		print OUT "\n$_";
	}#endif
	else {
		$_=~s/\s+//g;
		print OUT $_;
	}
}#endwhile
close IN;
close OUT;