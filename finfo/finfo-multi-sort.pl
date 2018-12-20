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
use lib $FindBin::RealBin;
use Fahrinfo_ubuntu;

### set variables
print_exit() unless ($#ARGV >= 0);

my $nrhaltestellen = join('',$ARGV[0]);
my $nroarg = 2 * $nrhaltestellen + 2;
print_exit() unless ($#ARGV == $nroarg);
my @haltestellen;
my @xnumm;
my $maxj = 500;
my @haltestellennr;
my @station;
my @temparr;
my @url;
my @command;
foreach my $xh (0..($nrhaltestellen - 1)){
 $haltestellen[$xh] = join('',$ARGV[($xh + 1)]);
 $xnumm[$xh] = join('',$ARGV[($xh + $nrhaltestellen + 2)]);
 @temparr = haltnummer($haltestellen[$xh],$xnumm[$xh]);
 push(@station,$temparr[0]);
 push(@haltestellennr,Fahrinfo::get_number($temparr[1]));
}
my $typus = join('',$ARGV[($nrhaltestellen + 1)]);
print_exit() if ($typus ne 'dep' and $typus ne 'arr');

my $filtre = join('',$ARGV[($nrhaltestellen + 2)]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo::calc_filter($filtre);

my $num;

foreach my $xh (0..$#haltestellennr) {
 $url[$xh] = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellennr[$xh]&boardType=$typus&maxJourneys=$maxj&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&";
 $command[$xh] = "elinks -dump -dump-width 200 \"$url[$xh]\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g;s/è/e/g;s/ê/e/g;s/à/a/g;s/á/a/g;s/â/a/g\'";
}


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
  my @pretext; 
  my @prenr;
  # open a pipe to the acpi command and read the battery value
  # and a few other parameters
  foreach my $xh (0..$#haltestellennr) {
   open(ACPI, "$command[$xh] |") || die "can't open pipe $xh!";
   Fahrinfo::read_str_pn(\*ACPI,$xh,$station[$xh],\@prenr,\@pretext);
   close ACPI;
  }
  Fahrinfo::st_entr_n($#haltestellennr,\@prenr,\@pretext,\@text);
  Fahrinfo::pmax_sgl(\@text);
  Fahrinfo::add_lbr(\@text);
  my $nrtext = $#text;
  if($nrtext > -1){
    print "@text" if ($debug);
  # set the labels and battery indicator
    $label->configure(-text=>"@enctext");
  }
}

sub haltnummer {
  my ($halt, $xnum) = @_;
  open(DATA, "$FindBin::RealBin/../data/fahrinfo-elinks2.dat") || die "can't open 'fahrinfo-elinks2.dat'";
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
    print "$haltarr[$numma]\n";
    @testarr = ($haltarr[$numma],$haltnrarr[$numma]);
    return @testarr;
  }elsif($numma > 0 and defined $xnum){
    print "$xnum : $haltarr[$xnum]\n";
    @testarr = ($haltarr[$xnum],$haltnrarr[$xnum]);
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
  print "Gebrauch: fahrinfo-multi-sort.pl <#Haltestellen> \"Haltestelle1\" ... <dep/arr> <code> <Zeile1> ... \n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}
