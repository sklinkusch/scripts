#!/usr/bin/perl 

use strict;                            # all variables have to be declared before use
use warnings;                          # print all errors and warnings to the shell
use utf8;                              # use unicode table for font encryption
use Encode;                            # packages to reencode text
use open ':encoding(utf8)';                                                      
use open ':std';                                                                 
use FindBin;                           # package to find an extra file for modules
use lib $FindBin::RealBin;                 # absolute path where extra file is found
use Fahrinfo_ubuntu;                   # name of extra file (without .pm ending)

### set variables
print_exit() unless ($#ARGV == 7);

my $ahaltestelle = join('',$ARGV[0]);
my $bhaltestelle = join('',$ARGV[4]);
my $axnumm = join('',$ARGV[3]);
my $bxnumm = join('',$ARGV[7]);
my @ahaltestellenarr = haltnummer($ahaltestelle,$axnumm);
my @bhaltestellenarr = haltnummer($bhaltestelle,$bxnumm);
my $ahaltestellennum = $ahaltestellenarr[1];
my $bhaltestellennum = $bhaltestellenarr[1];
my $ahaltestellenr = Fahrinfo_ubuntu::get_number($ahaltestellennum);
my $bhaltestellenr = Fahrinfo_ubuntu::get_number($bhaltestellennum);
my $astation = $ahaltestellenarr[0];
my $bstation = $bhaltestellenarr[0];

my $atypus = join('',$ARGV[1]);
print_exit() if ($atypus ne 'dep' and $atypus ne 'arr');
my $btypus = join('',$ARGV[5]);
print_exit() if ($btypus ne 'dep' and $btypus ne 'arr');

my $afiltre = join('',$ARGV[2]);
$afiltre = 127 if ($afiltre < 1 or $afiltre > 127);
my $afilter = Fahrinfo_ubuntu::calc_filter($afiltre);
my $bfiltre = join('',$ARGV[6]);
$bfiltre = 127 if ($bfiltre < 1 or $bfiltre > 127);
my $bfilter = Fahrinfo_ubuntu::calc_filter($bfiltre);

my $num;

my $aurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$ahaltestellenr&boardType=$atypus&maxJourneys=200&selectDate=today&productsFilter=$afilter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $burl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$bhaltestellenr&boardType=$btypus&maxJourneys=200&selectDate=today&productsFilter=$bfilter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
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
  my $num = shift;
  my @text; 
  my @atext;  
  my @aptext; 
  my @btext; 
  my @bptext; 
  open (ACPI, "$acommand |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_str(\*ACPI,$astation,\@aptext);
  close ACPI;
  Fahrinfo_ubuntu::make_fstr(\@aptext,\@atext);

  open (ACPJ, "$bcommand |") || die "can't open pipe!";
  Fahrinfo_ubuntu::read_str(\*ACPJ,$bstation,\@bptext);
  close ACPJ;
  Fahrinfo_ubuntu::make_fstr(\@bptext,\@btext);

  Fahrinfo_ubuntu::pmax_dual(\@atext,\@btext);
  Fahrinfo_ubuntu::cmp_fstr_dual(\@atext,\@btext,\@text);
  my $nrtext = $#text;
  if($nrtext > -1){
    print "@text";
  }
}

sub haltnummer {
  my $halt = $_[0];
  my $numm = $_[1];
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
    @testarr = ($haltarr[$numma],$haltnrarr[$numma]);
    return @testarr;
  }elsif($numma > 0 and defined $numm){
    @testarr = ($haltarr[$numm],$haltnrarr[$numm]);
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
  print "Gebrauch: finfo-fstrict-duo.pl \"Haltestelle1\" <dep/arr> <code> <Zeile> \"Haltestelle2\" <dep/arr> <code> <Zeile>\n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}
