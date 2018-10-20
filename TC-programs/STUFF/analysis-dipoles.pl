#!/usr/bin/perl -w

use strict;
use warnings;

my $file = join('',$ARGV[0]);
my $nros;
my $cvar;
my @energies;
my @ion;
my @dx;
my @dy;
my @dz;
my @dipol;
my $count;
my $icount;
my $dcount;
my $dumen;
my @dumdip;

$cvar = 0;
$count = 0;
$icount = 0;
$dcount = 0;
open my $A, $file or die "Cannot open $file: $!";
while(my $line = <$A>){
    if($line =~ /Nbasis/ and $cvar == 0){
	$nros = parsebasis($line);
	$cvar = 1;
    }
    if($line =~ /([\d]{1,3}[ \t]+[\d]{1,2}[\.]{0,1}[\d]{0,6})/ and $count < $nros){
	$dumen = parseens($line);
	$energies[$count] = $dumen;
	$count++;
	if($count == $nros){
	    $count++;
	    next;
	}
    }
    if($line =~ /([\d]{1,3}[ \t]+[\d]{1}[\.]{0,1}[\d]{0,6})/ and $count >= $nros and $icount < $nros){
	$dumen = parseens($line);
	$ion[$icount] = $dumen;
	$icount++;
    }
    if($line =~ /(0[ \t]+[\d]{1,3}[ \t]+[\-\d\.e]+[ \t]+[\-\d\.e]+[ \t]+[\-\d\.e]+)/ and $dcount < $nros){
	@dumdip = parsedip($line);
	$dx[$dcount] = $dumdip[0];
	$dy[$dcount] = $dumdip[1];
	$dz[$dcount] = $dumdip[2];
	$dipol[$dcount] = sqrt(($dx[$dcount]*$dx[$dcount])+($dy[$dcount]*$dy[$dcount])+($dz[$dcount]*$dz[$dcount]));
	$dcount++;
    }
}

foreach my $x(0..($nros-1)){
    printf("%.3u   %9.6f   %.6f   % .6f   % .6f    % .6f   %.6f\n",$x,$energies[$x],$dipol[$x],$dx[$x],$dy[$x],$dz[$x],$ion[$x]);
#    print "$x   $energies[$x]   $dipol[$x]   $dx[$x]   $dy[$x]   $dz[$x]   $ion[$x]\n";
}

sub parsebasis {
    my $str = shift;
    $str =~ /Nbasis:[ \t]+([\d]+)/;
    my $num = $1;
    return $num;
}

sub parseens {
    my $str = shift;
    $str =~ /([\d]{1,3})[ \t]+([\d]{1,2}[\.]{0,1}[\d]{0,6})/;
    my $num = $2;
    return $num;
}

sub parsedip {
    my $str = shift;
    $str =~ /(0)[ \t]+([\d]{1,3})[ \t]+([\-\d\.e]+)[ \t]+([\-\d\.e]+)[ \t]+([\-\d\.e]+)/;
    my @num;
    $num[0] = $3;
    $num[1] = $4;
    $num[2] = $5;
    return @num;
}

