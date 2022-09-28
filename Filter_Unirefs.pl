#!/usr/bin/perl

use Data::Dumper;
use strict;

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot open the metadata PID file\n";
my $hash={};

while (chomp (my $line = <IN>))
{
	my @meta = split(",", $line);
	my $sid = $meta[0];
	my $pid = $meta[1];
	$hash ->{$sid}=$pid;
}
print Dumper $hash;

my $input1 = <$ARGV[1]>;
open (IN1, $input1) || die "Cannot open the file Count matrix file\n";
chomp (my $line1 = <IN1>)
my @header = split(",", $line1);
print @header;

#while (chomp (my $line2 = <IN1>))
#{
#	my @uniref = split(",", $line2);
#	for my $i (@uniref)
#	{
#		if ($i > 0)
#		{
#			if (exists($hash->{$i}))
#			{
#				$hash->{$i} = 1;
#			}
#			else{
#				next;
#			}
#	
#		}
#	}	
#}
