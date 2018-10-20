#!/usr/bin/perl -w

use strict;
use warnings;

if($#ARGV != 1){
    print "Usage: ./get_spectraldata <GAMESS LOG> <output file>\n";
    exit;
}
my $infile = join('',$ARGV[0]);
if(-e $infile){
    if(-r $infile){
    }else{
	print "GAMESS logfile is not readable\n";
	exit;
    }
}else{
    print "GAMESS logfile is not existing\n";
    exit;
my $outfile = join('',$ARGV[1]);
if(-e outfile){
    if(-w "Output file is not writable\n";
    }
}
my $gsenergy;
my @exenergy;
my @osz;
my $gscount = 0;
my $excount = 0;

open my $log, $infile or die "Could not open $infile: $!";
open OUF, "> $outfile" or die "Could not open $outfile: $!";

while (my $line = <$log>) {
    if ($line =~ /TOTAL ENERGY =/) {
	$gsenergy = parseFloat($line);
	$gscount++;
    }
    if($line =~ /[\d]+[ ]{2}A[ ]+[-]*[\d]+[\.][\d]+[ ]+[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[\d]+[\.][\d]+/) {
	my $buf = parseLine($line);
	push @exenergy, $buf;
	my $buff = parseOsz($line);
	push @osz, $buff;
	$excount++;
    }
}
close($log);
if($gscount <= 0){
    print "No ground state energy found\n";
    exit;
}
if($excount <= 0){
    print "No excited states found\n";
    exit;
}

foreach my $x (@exenergy) {
    $x = $x - $gsenergy;
}
unshift(@exenergy,0.);
unshift(@osz,0.);
foreach my $x (1..$#exenergy) {
    print OUF "$exenergy[$x]    $osz[$x]\n";
}

sub parseFloat {
    my $str = shift;
    $str =~ /([-]*[\d\.]+[ ]*[-]*[\d\.]+[ ]*[-]*[\d\.]+)/;
    return $1;
}

sub parseLine {
    my $str = shift;
    $str =~ /[\d]+[ ]{2}A[ ]+([-]*[\d.]+)[ ]+[\d.]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+([\d]+[\.][\d]+)/;
    return $1;
}

sub parseOsz {
    my $str = shift;
    $str =~ /[\d]+[ ]{2}A[ ]+([-]*[\d.]+)[ ]+[\d.]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+([\d]+[\.][\d]+)/;
    return $2;
}

