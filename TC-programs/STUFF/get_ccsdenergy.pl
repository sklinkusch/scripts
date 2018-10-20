#!/usr/bin/perl

use strict;
use warnings;

my $file = join('', $ARGV[0]);
my $debtoau = 2.5417464;
my $gsenergy = 0.;
my @exhenergy;
my @exevenergy;
my @totalenergy;
my $stateno = 1;

open my $log, $file or die "Count not open $file: $!";

while (my $line = <$log>) {
    if ($line =~ /COUPLED-CLUSTER ENERGY/) {
	$gsenergy = parseFloat($line);
	$exhenergy[0] = 0.000;
	$exevenergy[0] = 0.000;
	$totalenergy[0] = $gsenergy;
    }
    if($line =~ /A[ \t]+(\d\.[\d]{1,8})[ \t]+([\d]{1,3}\.[\d]{1,4})[ \t]+([-]{0,1}[\d]{1,4}\.[\d]{1,8})[ \t]+CONVERGED/) {
	my @buf = split(/[ ]+/, parseLine($line));
	$exhenergy[$stateno] = $buf[0];
	$exevenergy[$stateno] = $buf[1];
	$totalenergy[$stateno] = $buf[2];
        $stateno++;
    }
    if($line =~ /A[ \t]+(\d\.[\d]{1,8})[ \t]+([\d]{1,3}\.[\d]{1,4})[ \t]+([-]{0,1}[\d]{1,4}\.[\d]{1,8}[ \t])+NOT CONVERGED/) {
	$stateno++;
    }
}
close $log;

foreach my $x (0..$#exhenergy){
    if(defined $totalenergy[$x]){
	printf "%2d  % 13.8f  %11.8f  %6.3f\n", $x,  $totalenergy[$x], $exhenergy[$x], $exevenergy[$x];
    }else{
	print "Total energy for $x not defined\n";
    }
}

sub parseFloat {
    my $str = shift;
    $str =~ /([-]{0,1}[\d\.]+)/;
    return $1;
}

sub parseLine {
    my $str = shift;
    $str =~ /([\d]{1,2}[\.]{1}[\d]{1,8})[ \t]+([\d]{1,2}[\.]{1}[\d]{1,3})[ \t]+([-]{0,1}[\d]{1,4}[\.]{1}[\d]{1,8})/;
    my $num;
    $num = "$1 $2 $3";
    return $num;
}
