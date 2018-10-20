#!/usr/bin/perl -w

use strict;
use warnings;

my $file = join('',$ARGV[0]);
my @string;
my %excitation;

open my $A, $file or die "Cannot open $file: $!\n";
while(my $line = <$A>){
    if($line =~ /([\d]{1,3}[ \t]+[\d]{1,3}[ \t]+[\-]*[\d]{1}[\.]{1}[\d]{6}[ \t]+[\-]*[\d]{1}[\.]{1}[\d]{6})/){
	@string = parseline($line);
	$excitation{$string[0]} = $string[1];
    }
}
close $A;

foreach my $key (sort { $excitation{$b} <=> $excitation{$a} } keys %excitation){
    printf("%12s : %14.12f\n", $key, $excitation{$key});
}

sub parseline {
    my $str = shift;
    $str =~ /([\d]{1,3})[ \t]+([\d]{1,3})[ \t]+[\-]*([\d]{1}[\.]{1}[\d]{6})[ \t]+[\-]*[\d]{1}[\.]{1}[\d]{6}/;
    my $occ = $1;
    my $virt = $2;
    my $exc = "$occ -> $virt";
    my $num = $3;
    my @numm;
    $numm[0] = $exc;
    $numm[1] = $num;
    return @numm;
}
