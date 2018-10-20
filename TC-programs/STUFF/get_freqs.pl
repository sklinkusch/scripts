#!/usr/bin/perl -w

use strict;
use warnings;

if($#ARGV != 0){
    print "Usage: ./get_freqs.pl <gamess log>\n";
    exit;
}
my $file = join('',$ARGV[0]);
if(-e $file){
    if(-r $file){
    }else{
	print "File could not be read\n";
	exit;
    }
}else{
    print "File does not exist\n";
    exit;
}
my @freq;
my @ints;
my $dumpf;
my $transrot;

open my $A, $file or die "Could not open $file: $!";
while (my $line = <$A>){
    if($line =~ /MODES 1 TO [\d]+ ARE TAKEN AS ROTATIONS AND TRANSLATIONS./){
	$transrot = parseint($line);
    }
    if($line =~ /FREQUENCY:[ \t]+[\d\.]+[ \t]*[\d\.]*[ \t]*[\d\.]*[ \t]*[\d\.]*[ \t]*[\d\.]*/){
	my @dumfreq = parseFreq($line);
	my $maxfreq = $dumfreq[5] - 1;
	foreach my $x (0..$maxfreq){
	    push (@freq, $dumfreq[$x]);
	}
    }
    if($line =~ /IR INTENSITY:[ \t]+[\d\.]+[ \t]*[\d\.]*[ \t]*[\d\.]*[ \t]*[\d\.]*[ \t]*[\d\.]*/){
	my @dumints = parseInts($line);
	my $maxints = $dumints[5] - 1;
	foreach my $x (0..$maxints){
	    push (@ints, $dumints[$x]);
	}
    }
}
close $A;
if($#freq != $#ints){
    print "Error while parsing: number of freqs and ints not identical\n";
    exit;
}
for(my $x = 0; $x < $transrot; $x++){
    shift @ints;
    shift @freq;
}

foreach my $x (0..$#freq){
    $freq[$x] /= 219474.83;
    printf("%12.10f   %7.5f\n",$freq[$x], $ints[$x]);
}

sub parseFreq {
    my $str = shift;
    $str =~ /FREQUENCY:[ \t]+([\d\.]+)[ \t]*([\d\.]*)[ \t]*([\d\.]*)[ \t]*([\d\.]*)[ \t]*([\d\.]*)/; 
    my @numm;
    $numm[0] = $1;
    $numm[5] = 1;
    if(defined $2){
        $dumpf = $2;
	if($dumpf =~ /[\d\.]+/){
	    $numm[1] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIN;
	}
    }else{
	goto FIN;
    }
    if(defined $3){
	$dumpf = $3;
	if($dumpf =~ /[\d\.]+/){
	    $numm[2] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIN;
	}
    }else{
	goto FIN;
    }
    if(defined $4){
	$dumpf = $4;
	if($dumpf =~ /[\d\.]+/){
	    $numm[3] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIN;
	}
    }else{
	goto FIN;
    }
    if(defined $5){
	$dumpf = $5;
	if($dumpf =~ /[\d\.]+/){
	    $numm[4] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIN;
	}
    }else{
	goto FIN;
    }
    FIN:
    return @numm;
}

sub parseInts {
    my $str = shift;
    $str =~ /IR INTENSITY:[ \t]+([\d\.]+)[ \t]*([\d\.]*)[ \t]*([\d\.]*)[ \t]*([\d\.]*)[ \t]*([\d\.]*)/;
    my @numm;
    $numm[0] = $1;
    $numm[5] = 1;
    if(defined $2){
	$dumpf = $2;
	if($dumpf =~ /[\d\.]+/){
	    $numm[1] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIA;
	}
    }else{
	goto FIA;
    }
    if(defined $3){
	$dumpf = $3;
	if($dumpf =~ /[\d\.]+/){
	    $numm[2] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIA;
	}
    }else{
	goto FIA;
    }
    if(defined $4){
	$dumpf = $4;
	if($dumpf =~ /[\d\.]+/){
	    $numm[3] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIA;
	}
    }else{
	goto FIA;
    }
    if(defined $5){
	$dumpf = $5;
	if($dumpf =~ /[\d\.]+/){
	    $numm[4] = $dumpf;
	    $numm[5] += 1;
	}else{
	    goto FIA;
	}
    }else{
	goto FIA;
    }
    FIA:
    return @numm;
}

sub parseint {
    my $str = shift;
    $str =~ /MODES 1 TO ([\d]+) ARE TAKEN AS ROTATIONS AND TRANSLATIONS./;
    my $num = $1;
    return $num;
}
