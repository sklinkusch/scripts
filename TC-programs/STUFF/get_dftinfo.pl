#!/usr/bin/perl

use strict;
use warnings;

if($#ARGV != 0){
    print "Usage: ./get_dftinfo.pl <GAMESS LOG> > output.dat\n";
    exit;
}
my $file = join('', $ARGV[0]);
if(-e $file){
    if(-r $file){
    }else{
	print "Input file is not readable\n";
	exit;
    }
}else{
    print "Input file does not exist\n";
    exit;
}
my $debtoau = 2.5417464;
my $gsenergy = 0.;
my @exenergy;
my $stateno = 0;
my @dipx;
my @dipy;
my @dipz;
my @buf;
my @states;
my @frommo;
my @tomo;
my @uprate;
my @downrate;
my $oldstateno = 0;
my $statecount = 0;
my @nostates;
my $dumvar;
my $y = 0;
my @substatus;
my @subfrommo;
my @subtomo;
my @subuprate;
my @subdnrate;
my $gscount = 0;
my $excount = 0;

open my $log, $file or die "Could not open $file: $!";

while (my $line = <$log>) {
    if ($line =~ /TOTAL ENERGY =/) {
	$gsenergy = parseFloat($line);
	$gscount++;
    }
    if($line =~ /[\d]+[ ]{2}A[ ]+[-]*[\d]+[\.][\d]+[ ]+[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[\d]+[\.][\d]+/) {
	my @buf = split(/[ ]+/, parseLine($line));
	push @exenergy, $buf[2];
    }
    if($line =~ /STATE #.*ENERGY/ or $line =~ /[0-9]{1,3}[ ]{4,7}[0-9]{1,3}[ ]{3,4}[\-]{0,1}[0-9]{1}\.[0-9]{6}[ ]{6,7}[\-]{0,1}[0-9]{1}\.[0-9]{6}/) {
	if($line =~ /STATE #.*ENERGY/){
	    $stateno = parseInt($line);
	    $excount++;
	}
	if($line =~ /[0-9]{1,3}[ ]{4,7}[0-9]{1,3}[ ]{3,4}[\-]{0,1}[0-9]{1}\.[0-9]{6}[ ]{6,7}[\-]{0,1}[0-9]{1}\.[0-9]{6}/){
	    my @buf = split(/[ ]+/, parseEx($line));
	    if($stateno != 0) {
		push @states, $stateno;
		push @frommo, $buf[0];
		push @tomo, $buf[1];
		push @uprate, $buf[2];
		push @downrate, $buf[3];
		if($oldstateno == $stateno){
		    $statecount++;
		}elsif($oldstateno == $stateno - 1){
		    push(@nostates, $statecount);
		    $statecount = 1;
		}else{
		    my $pristate = $oldstateno + 1;
		    print "No eigenvector elements found for eigenstate $pristate\n";
		    exit;
		}
		$oldstateno = $stateno;
	    }
	}
    }
}
if($gscount == 0){
    print "No ground state energy found in input file\n";
    exit;
}
if($excount == 0){
    print "No excitation energies found in input file\n";
    exit;
}
push(@nostates, $statecount);
my $count = 1;
foreach my $x (@exenergy) {
    $x = $x - $gsenergy;
    $count++;
}
unshift(@exenergy,0.);
if($states[$#states] != $#exenergy){
    print "Nr of states not equal\n";
    exit;
}    
my $nros = @exenergy;

print $gsenergy, "\n";
print $nros, "\n";
foreach my $x (@exenergy) {
    print $x, "\n";
}
for(my $x = 0; $x < $nros; $x++){
    print $nostates[0], "\n";
    $dumvar = $nostates[0];
    shift(@nostates);
    $y = 0;
    while($y < $dumvar){
	print $states[0], " ", $frommo[0], " ", $tomo[0], " ", $uprate[0], " ", $downrate[0], "\n";
	$y++;
	shift(@states);
	shift(@frommo);
	shift(@tomo);
	shift(@uprate);
	shift(@downrate);
    }
}

close $log;

sub parseFloat {
    my $str = shift;
    $str =~ /([-]*[\d\.]+[ ]*[-]*[\d\.]+[ ]*[-]*[\d\.]+)/;
    return $1;
}

sub parseInt {
    my $str = shift;
    $str =~ /([\d]+)/;
    return $1;
}

sub parseEx {
    my $str = shift;
    $str =~ /([0-9]{1,3}[ ]{4,7}[0-9]{1,3}[ ]{3,4}[\-]{0,1}[0-9]{1}\.[0-9]{6}[ ]{6,7}[\-]{0,1}[0-9]{1}\.[0-9]{6})/;
    return $1;
}

sub parseLine {
    my $str = shift;
    $str =~ /([\d]+[ ]{2}A[ ]+[-]*[\d]+[\.][\d]+[ ]+[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[-]*[\d]+[\.][\d]+[ ]+[\d]+[\.][\d]+)/;
    return $1;
}
