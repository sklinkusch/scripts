#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

my $file = join("", $ARGV[0]);
my @dataarray = ();
my $lineno = 0;
my @sourcearray = ("Firma", "Anrede", "Titel", "Vorname", "Nachname", "Strasse", "Hausnummer", "Postfach", "PLZ", "Ort");
open(CSVIN, $file) || die "Cannot open $file!\n";
while (my $line = <CSVIN>) {
  if($lineno == 0 || $line =~ /Packstation/) {
    $lineno++;
    next;
  }
  my %obj;
  my @arr = split(/,/, $line);
  foreach my $ioarr (0..$#arr) {
    if (length($arr[$ioarr]) > 0) {
      $obj{$sourcearray[$ioarr]} = $arr[$ioarr];
    }
  }
  push(@dataarray, \%obj);
  $lineno++;
}
close(CSVIN);

my @only_germany = ();
my @not_germany = ();
foreach my $zz (0..$#dataarray) {
  push(@only_germany, $dataarray[$zz]) if($dataarray[$zz]->{PLZ} =~ /^[0-9]{5}$/);
}

my @final;
foreach my $iodarr (0..$#only_germany) {
  my $streetsrc = $only_germany[$iodarr]->{Strasse};
  my ($street, $number);
  my ($firma, $anrede, $titel, $vorname, $nachname, $postfach, $plz, $ort);
  if(defined $only_germany[$iodarr]->{Firma}) {
      $firma = $only_germany[$iodarr]->{Firma};
    } else {
      $firma = "";
    }
    if(defined $only_germany[$iodarr]->{Anrede}) {
      $anrede = $only_germany[$iodarr]->{Anrede};
    } else {
      $anrede = "";
    }
    if(defined $only_germany[$iodarr]->{Titel}) {
      $titel = $only_germany[$iodarr]->{Titel};
    } else {
      $titel = "";
    }
    if(defined $only_germany[$iodarr]->{Vorname}) {
      $vorname = $only_germany[$iodarr]->{Vorname};
    } else {
      $vorname = "";
    }
    if(defined $only_germany[$iodarr]->{Nachname}) {
      $nachname = $only_germany[$iodarr]->{Nachname};
    } else {
      $nachname = "";
    }
    if(defined $only_germany[$iodarr]->{Postfach}) {
      $postfach = $only_germany[$iodarr]->{Postfach};
    } else {
      $postfach = "";
    }
    if(defined $only_germany[$iodarr]->{PLZ}) {
      $plz = $only_germany[$iodarr]->{PLZ};
    } else {
      $plz = "";
    }
    if(defined $only_germany[$iodarr]->{Ort}) {
      $ort = $only_germany[$iodarr]->{Ort};
    } else {
      $ort = "";
    }
  if($streetsrc =~ /^[\s]*([A-ZÄÖÜa-zäöüß\-\.\'\/ ]+)[ ]?([0-9\/\-\+ ]+[a-zA-Z]?)[\s]*$/){
    ($street, $number) = separate_address($streetsrc);
    printf "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s", $firma, $anrede, $titel, $vorname, $nachname, $street, $number, $postfach, $plz, $ort;
  } else {
    if(defined $only_germany[$iodarr]->{Strasse}) {
      $street = $only_germany[$iodarr]->{Strasse};
    } else {
      $street = "";
    }
    if(defined $only_germany[$iodarr]->{Hausnummer}) {
      $number = $only_germany[$iodarr]->{Hausnummer};
    } else {
      $number = "";
    }
    printf "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s", $firma, $anrede, $titel, $vorname, $nachname, $street, $number, $postfach, $plz, $ort;
  }
}

sub separate_address {
  my $str = shift;
  $str =~ /^[\s]*([A-Za-zäöüßÄÖÜ\-\.\'\/\s]+)[\s+]?([0-9\/\-\+\s]+[a-zA-Z]?)[\s]*$/;
  my $street = $1;
  my $number = $2;
  return ($street, $number);
}
