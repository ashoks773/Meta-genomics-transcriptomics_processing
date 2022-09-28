#!/usr/bin/perl

use Data::Dumper;
use strict;

use vars qw($opt_I $opt_C);
use Getopt::Std;
getopts("I:C:");

my $i=0;
my $ln = '';
my %hash=();

@ARGV=($opt_I);
while($ln = <>) {
	$i+=1;
	chomp($ln);
	$hash{$i} = $ln;
}

#print Dumper $hash;

@ARGV=($opt_C);
my @templine = ();
my %presence = ();
my $npatients = 0;
my $nsamples = 0;
my $ncounts = 0;

while($ln = <>) {
@templine=split(",",$ln);
my $records = scalar(@templine);
for(my $i=0; $i<$records; $i++) {
	if ($templine[$i]>0) {				### we can change "0" for a pre-defined threshold, to retrieve genes present in patients with at least N counts in any sample from that patient
		$ncounts +=$templine[$i];		### also, adding over $templine gives the total counts for each gene
		$nsamples +=1;				### also, easy to retrieve sample presence, along with patient presence

		if (defined $presence{$hash{$i+1}}) {	
			next;
		}
		else {
			$presence{$hash{$i+1}}=1;
			$npatients +=1;
		}
	}
}
print STDOUT $ncounts."\t".$nsamples."\t".$npatients."\n";
$npatients=0;
$nsamples=0;
$ncounts=0;
%presence = ();
}

