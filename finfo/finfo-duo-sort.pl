#!/usr/bin/perl 

#
# A Perl/Tk wrapper to graphically show the arrival and departure table for a certain station/tramway stop/bus stop
#
# It calls 'elinks', takes the content of the page, and writes it to # the graphical window. It is further updated every 10 seconds.
#
# Written by SK
#

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
print_exit() unless ($#ARGV == 5);

my $ahaltestelle = join('',$ARGV[0]);
my $bhaltestelle = join('',$ARGV[1]);
my $axnumm = join('',$ARGV[4]);
my $bxnumm = join('',$ARGV[5]);
my $xnumm = $axnumm;
my @ahaltestellenarr = haltnummer($ahaltestelle);
$xnumm = $bxnumm;
my @bhaltestellenarr = haltnummer($bhaltestelle);
my $ahaltestellenum = $ahaltestellenarr[1];
my $bhaltestellenum = $bhaltestellenarr[1];
my $ahaltestellenr = Fahrinfo_ubuntu::get_number($ahaltestellenum);
my $bhaltestellenr = Fahrinfo_ubuntu::get_number($bhaltestellenum);
my $astation = $ahaltestellenarr[0];
my $bstation = $bhaltestellenarr[0];

my $typus = join('',$ARGV[2]);
print_exit() if ($typus ne 'dep' and $typus ne 'arr');

my $filtre = join('',$ARGV[3]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo_ubuntu::calc_filter($filtre);

my $num;

my $aurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$ahaltestellenr&boardType=$typus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $burl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$bhaltestellenr&boardType=$typus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $acommand = "elinks -dump -dump-width 200 \"$aurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $bcommand = "elinks -dump -dump-width 200 \"$burl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

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
  my ($num) = @_;
  my (@text);
  my @apretext; 
  my @bpretext; 
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  open (ACPI, "$acommand |") || die "can't open pipe A!";
  Fahrinfo_ubuntu::read_str_p(\*ACPI,$astation,\@apretext);
  close ACPI;
  open (ACPJ, "$bcommand |") || die "can't open pipe B!";
  Fahrinfo_ubuntu::read_str_p(\*ACPJ,$bstation,\@bpretext);
  close ACPJ;
  Fahrinfo_ubuntu::st_entr(\@apretext,\@bpretext,\@text);
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
  print "Gebrauch: finfo-duo-sort.pl \"Haltestelle1\" \"Haltestelle2\" <dep/arr> <code> <Zeile1> <Zeile2>\n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}
