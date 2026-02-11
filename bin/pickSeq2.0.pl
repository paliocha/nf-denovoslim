#!/usr/bin/perl
#Last modified 22.8.2011
#2.8.2011 version 0.1
#22.8.2011 version 0.2 also header search without additional information
use strict;
use Getopt::Long;
my $VERSION=0.2;

my $helpMsg="
********************************************************************************
pickSeq.pl

Extracts reads with certain identifiers from fastq or fasta file.

Input:
-i string   Fasta or fastq illumina1.3+ file to be searched. Required.
-q string   Sequence headers to be queried, one header per row. Optional.
-s string   One header. Optional.
            Example inputs: 
            myheader
            \"myheader some additional info text\"
			\@fastqheader
            Use parameters -q and/or -s to define the query header.
-o string   output fasta file where found sequences are printed. Required.
-t          trim fasta file headers so additional information isn't searched.

            Examples:
            query          fastafile header   match without -t  match with -t
            -----          ----------------   ----------------  -------------
            myheader       >myheader          YES               YES
            myheader       >myheader info     NO                YES
            myheader info  >myheader          NO                NO
            myheader	
		
-h          help message

-v          print all that is not found

#version 2.0: -x removed, -v added
#version 0.2: whitespace characters in query id:s and fasta sequence header 
#             ends are trimmed
********************************************************************************
";
my $infile;
my $query;
my $outfile;
my $help;
my $queryids=0;
my $trimheader;
my $found=0;
my $seqheader;
my $inversesearch=0;
GetOptions("i:s"=> \$infile,"q:s"=> \$query,"s:s"=> \$seqheader, "o:s"=> \$outfile,"t"=> \$trimheader,"h"=> \$help,"v"=>\$inversesearch); 
if($help==1) { print $helpMsg; exit; }	

	print"Reading query id:s...\n";
	my %ids=();	
	if($query ne "") {
	open(QUERY,$query) or die("Give queried sequence id:s as param -q!");
	while(<QUERY>) {
		chomp;
		$_=~s/\s+$//;
		print "$_\n";
		$ids{$_}=1; 
		if(length($_)>0) { $queryids++; }
	}#endwhile
}#endif
if($seqheader ne "") { $ids{$seqheader}=1; $queryids=1; print "$seqheader\n";  }
close QUERY;

print("Infile: $infile, outfile: $outfile.\n");
my $len=keys(%ids);

#check whether the infile is in fasta or fastq format
open(INFILE,$infile) || die("can't open infile $infile");
my $line=<INFILE>;
my $fasta=-1;
if($line=~m/^>/)  {$fasta=1;}
elsif($line=~m/^@/)  {$fasta=0;}
else {die("unknown file format in $infile");}
close INFILE;

#search string for a fasta header
my $searchstring='^>.+';
if ($trimheader==1) { $searchstring='^>\S+'; }  

my $id;
my $id_trimmed;
open(INFILE, "$infile") or die("Can't open input file $infile. Give input file as param -i.");
open(OUTFILE,">$outfile") or die("Can't open output file $outfile Give output file as param -o.");

my $continuewriting='false';

#search fastq format
if($fasta==0) {
	while($line=<INFILE>) {
		$line=~s/\s*$//;
		$id=$line;
		
		if(($inversesearch==0 & exists($ids{$id})) || ($inversesearch==1 & !exists($ids{$id}))) {
			print OUTFILE $line."\n";
			my $line= <INFILE>;
			print OUTFILE $line;
			my $line= <INFILE>;
			print OUTFILE $line;
			my $line= <INFILE>;
			print OUTFILE $line;
			$found++;
		}
		else {<INFILE>;<INFILE>;<INFILE>;} 
	}
}
if($fasta==1) {
	while($line=<INFILE>) {
			
		if($line=~m/($searchstring)/) { #if header line then 
			$id=$1;
			$id=~s/^>//; 
			$id_trimmed=$id;
			$id_trimmed=~s/\s+$//;
			
			my $item;
		
			if(($inversesearch==0 & exists($ids{$id})) || ($inversesearch==1 & !exists($ids{$id}))) { 

				$continuewriting='true';
				$found++;
				print OUTFILE "$line";
			}
			else { $continuewriting='false';}	

		}#endif
			
		else { #if non header line
			if($continuewriting eq 'true') {  print OUTFILE "$line"; }
		
		}#endelse (if non'header line
	}#endwhile
}
close(OUTFILE);
close(INFILE);
print"Found $found sequences using $queryids query id:s.";