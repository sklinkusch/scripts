#!/usr/bin/perl -w

use strict;
use warnings;

my %vecs;
my %avecs;
my $num;
my @numm;
my $vec;
my $file = join('',$ARGV[0]);
my $state = join('',$ARGV[1]);
my $vecsdim;

open my $A, $file or die "Could not open $file: $!\n";
while(my $line = <$A>){
    if($line =~ /([\d]+[ \t]+[\d]+[ \t]+[\d\.\+\-edE]+)/){
	@numm = parseline($line);
	if($numm[0] == $state){
	    $vecs{$numm[1]} = $numm[2];
	}else{
	    next;
	}
    }
}
close $A;

$vecsdim = scalar keys %vecs;
if ($vecsdim != 0){
    while ((my $key, my $val) = each %vecs){
	if($val < 0){
	    $avecs{$key} = -1.*$val;
	}else{
	    $avecs{$key} = $val;
	}
    }

    foreach my $key (sort { $avecs{$b} <=> $avecs{$a} } keys %avecs){
	printf("%5u   %5u   % 32.30f\n",$state,$key,$vecs{$key});
    }
}

sub parseline {
    my $str = shift;
    my @num;
    $str =~ /([\d]+)[ \t]+([\d]+)[ \t]+([\d\.\+\-edE]+)/;
    $num[0] = $1;
    $num[1] = $2;
    $num[2] = $3;
    return @num;
}

