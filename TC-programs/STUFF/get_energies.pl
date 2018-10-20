#!/usr/bin/perl

use strict;
use warnings;

if($#ARGV != 1){
    print "Usage: ./get_energies.pl <GAMESS LOG> <CIS(D) LOG> > output.ens\n";
    exit;
}
my $cisfile = join('', $ARGV[0]);
if(-e $cisfile){
    if(-r $cisfile){
    }else{
	print "GAMESS logfile not readable\n";
	exit;
    }
}else{
    print "GAMESS logfile does not exist\n";
    exit;
}
my $cisdfile = join('', $ARGV[1]);
if(-e $cisdfile){
    if(-r $cisdfile){
    }else{
	print "CIS(D) logfile not readable\n";
	exit;
    }
}else{
    print "CIS(D) logfile does not exist\n";
    exit;
}
my $debtoau = 2.5417464;
my $count = 0;
my $rhfenergy;
my @cisenergy = ();
my @dcorrect;
my $rhfcount = 0;

open my $cis, $cisfile or die "Count not open $cisfile: $!";

while (my $line = <$cis>) {
    if ($line =~ /FINAL RHF ENERGY IS/) {
	$rhfenergy = parseFloat($line);
	$rhfcount++;
    }
    if ($line =~ /EXCITED STATE[0-9 ]*ENERGY=/) {
	push @cisenergy, parseEnergy($line);
	$count++;
    }
}
close $cis;
if($rhfcount == 0 or $rhfcount > 1 or $count == 0){
    print "Data not found in GAMESS logfile\n";
    exit;
}
foreach my $x (@cisenergy) {
    $x = $x - $rhfenergy;
}

open my $cisd,  $cisdfile or die "Could not open $cisdfile: $!";

$count = 0;
while (my $line = <$cisd>) {
	if($line =~ /Total corr/){
	    push @dcorrect, parseFloat($line);
	    $count++;
	}
}
close $cisd;
if($count == 0){
    print "Data not found in CIS(D) logfile\n";
    exit;
}
for (my $i = 0; $i < $count; $i++) {
    printf "%4.10f   %2.7f\n", $cisenergy[$i], $dcorrect[$i];
}

sub parseFloat {
    my $str = shift;
    $str =~ /([-]*[\d\.]+)/;
    return $1;
}

sub parseEnergy {
    my $str = shift;
    $str =~/(-[\d]+\.[\d]+)/;
    return $1;
}
