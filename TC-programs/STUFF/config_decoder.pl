#!/usr/bin/perl -w

use strict;
use warnings;

if($#ARGV != 2){
    print "Usage: ./config_decoder.pl <llim> <ulim> <nroe>\n";
    exit;
}
my $llim = join('',$ARGV[0]);
if(int_test($llim) == 0){
    print "$llim is not a positive integer!\n";
    exit;
}
my $ulim = join('',$ARGV[1]);
if(int_test($ulim) == 0){
    print "$ulim is not a positive integer!\n";
    exit;
}
my $nroe = join('',$ARGV[2]);
if(int_test($nroe) == 0){
    print "$nroe is not a positive integer!\n";
    exit;
}
if($nroe%2 != 0){
    print "$nroe is not even!\n";
    exit;
}
my $omo = $nroe/2 - $llim;
if($omo <= 0){
    print "No occupied orbitals in active space!\n";
    exit;
}
my $umo = $ulim - $nroe/2 + 1;
if($umo <= 0){
    print "No virtual orbitals in active space!\n";
    exit;
}
my $cissize = $omo * $umo + 1;
my $occ;
my $virt;

printf("%5u   %3u   %3u\n",0,0,0);

for(my $x = 1; $x < $cissize; $x++){
    $occ = ($x-1)/$umo+$llim;
    $virt = ($x-1)%$umo+$omo+$llim;
    printf("%5u   %3u   %3u\n",$x,$occ,$virt);
}

sub int_test {
    my $str = shift;
    my $num;
    if($str =~ /^\d+$/){
	$num = 1;
    }else{
	$num = 2;
    }
    return $num;
}

