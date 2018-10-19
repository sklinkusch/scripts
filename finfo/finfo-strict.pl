#!/usr/bin/perl 

use strict;
use warnings;
use utf8;
use Encode;
use open ':encoding(utf8)';
use open ':std';
use FindBin;
use lib $FindBin::Bin;
use Fahrinfo_ubuntu;

### set variables
print_exit() unless ($#ARGV == 2 or $#ARGV == 3);

my $haltestelle = join('',$ARGV[0]);
my $xnumm;
if($#ARGV == 3){
  $xnumm = join('',$ARGV[3]);
}
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellenum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo_ubuntu::get_number($haltestellenum);
my $station = $haltestellenarr[0];

my $typus = join('',$ARGV[1]);
print_exit() if ($typus ne 'dep' and $typus ne 'arr');

my $filtre = join('',$ARGV[2]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo_ubuntu::calc_filter($filtre);

my $num;

my $url = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$typus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $command = "elinks -dump -dump-width 200 \"$url\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

# call the power checker now
checkNet ($num);

############
# end main #
############

#
# subroutine to run the checking command, parse its output and 
# call the update subroutine
#
sub checkNet {
  my $num = shift;
  my (@text);
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$command |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_str(\*ACPI,$station,\@text);
  close ACPI;
  Fahrinfo_ubuntu::pmax_sgl(\@text);
  Fahrinfo_ubuntu::add_lbr(\@text);
  my $nrtext = $#text;
  if($nrtext > -1){
    print "@text";
  }
}

sub haltnummer {
  my $halt = shift;
  open(DATA, "$FindBin::Bin/../data/fahrinfo-elinks2.dat") || die "can't open 'fahrinfo-elinks2.dat'";
  my $numma = -1;
  my $str;
  my $haltu;
  my @haltarr;
  my @haltnrarr;
  while (my $line = <DATA>) {
    $haltu = decode_utf8($halt);
    $str = decode_utf8($line);
    if(index($line, $halt) != -1) {
      $numma++;
      $str =~ /(\"[a-zA-Z0-9\.\,\-\/\+\[\]\(\) ]+\")([\t])([0-9]+)/;
      push(@haltarr,$1);
      push(@haltnrarr,$3);
    }
  }
  my @testarr;
  if($numma == 0){
    @testarr = ($haltarr[$numma],$haltnrarr[$numma]);
    return @testarr;
  }elsif($numma > 0 and defined $xnumm){
    @testarr = ($haltarr[$xnumm],$haltnrarr[$xnumm]);
    return @testarr;
  }else{
    my $lnum = $numma + 1;
    print "$lnum Haltestellen gefunden: \n";
    foreach my $axnum (0..$#haltarr) {
      print "$axnum : $haltarr[$axnum]\n";
    }
    exit;
  }
}

sub print_exit {
  print "Gebrauch: finfo-strict.pl \"Haltestelle\" <dep/arr> <code> <Zeile>\n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}
