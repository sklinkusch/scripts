#!/usr/bin/perl -w

use strict;
use warnings;

my $filetot = join('',$ARGV[0]);
my $filea = join('',$ARGV[1]);
my $fileb = join('',$ARGV[2]);
my @string;
my %hashtot;
my %hasha;
my %hashb;
my %hashzu;

open my $tot, $filetot or die "Cannot open $filetot: $!";
while (my $line = <$tot>) {
    if($line =~ /[\d]{1,3}[ \t]+[\d][\.][\d]{12,20}[\-E\d]*/){
	@string = parseline($line);
	$hashtot{$string[0]} = $string[1];
    }
}
close $tot;

open my $strma, $filea or die "Cannot open $filea: $!";
while (my $line = <$strma>) {
    if ($line =~ /[\d]{1,3}[ \t]+[\d][\.][\d]{12,20}[\-E\d]*/){
	@string = parseline($line);
	$hasha{$string[0]} = $string[1];
    }
}
close $strma;

open my $strmb, $fileb or die "Cannot open $fileb: $!";
while (my $line = <$strmb>) {
    if ($line =~ /[\d]{1,3}[ \t]+[\d][\.][\d]{12,20}[\-E\d]*/){
	@string = parseline($line);
	$hashb{$string[0]} = $string[1];
    }
}
close $strmb;

for(my $x = 0; $x < 502; $x++){
    my $totval = $hashtot{$x};
    for(my $y = 0; $y < 251 ; $y++){
	if (defined $hasha{$y}){
	    my $aval = $hasha{$y};
	    if ($aval eq $totval) {
		$hashzu{$x} = "1,$y";
		delete $hasha{$y};
		last;
	    }
	}
    }
    for(my $y = 0; $y < 251; $y++){
	if(defined $hashb{$y}){
	    my $bval = $hashb{$y};
	    if($bval eq $totval) {
		$hashzu{$x} = "2,$y";
		delete $hashb{$y};
		last;
	    }
	}
    }
}

foreach my $x (0..501){
    if (defined $hashzu{$x}){
	print "$x : $hashzu{$x}\n";
    }else{
	print "$x : $hashtot{$x}\n";
    }
}

print "ISOMER 1:\n";
while((my $keya, my $vala) = each %hasha){
       print "1,$keya : $vala\n";
}
print "ISOMER 2:\n";
while((my $keyb, my $valb) = each %hashb){
       print "2,$keyb : $valb\n";
}

sub parseline {
    my $str = shift;
    $str =~ /([\d]{1,3})[ \t]+([\d]{1}[.]{1}[\d]{12,20}[E]{0,1}[-]{0,1}[\d]{3})/;
    my @numm;
    $numm[0] = $1;
    $numm[1] = $2;
    return @numm;
}

