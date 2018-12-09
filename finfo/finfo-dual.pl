#!/usr/bin/perl 

use strict;                         # all variables have to be declared before use
use warnings;                       # print all errors and warnings to the shell
use utf8;                           # use unicode table for font encryption
use Encode;                         # packages to reencode text
use open ':encoding(utf8)';
use open ':std';
use FindBin;                        # package to find an extra file for modules
use lib $FindBin::RealBin;              # absolute path where extra file is found
use Fahrinfo_ubuntu;                # name of extra file (without .pm ending), in brackets: loaded subroutines

### set variables
# exit code if nr of arguments is wrong
print_exit() unless ($#ARGV == 1 or $#ARGV == 2);

my $haltestelle = join('',$ARGV[0]);                               # input name of the station/bus stop
my $xnumm = join('',$ARGV[2]) if ($#ARGV == 2);                    # number defining the station/bus stop from a list of suitable stations/bus stops
my @haltestellenarr = haltnummer($haltestelle);                    # connect the name with the number
my $haltestellennum = $haltestellenarr[1];                         # preliminary number of the station/bus stop
my $haltestellenr = Fahrinfo_ubuntu::get_number($haltestellennum);        # correct number (insert '00' between first and second number)
my $station = $haltestellenarr[0];                                 # correct number of the station/bus stop

my $filtre = join('',$ARGV[1]);                                    # input number which products shall be used
$filtre = 127 if ($filtre < 1 or $filtre > 127);                   # correct a number that is out of range
my $filter = Fahrinfo_ubuntu::calc_filter($filtre);                # translate number to a code used in the URL

my $num;                                                           # dummy variable
my $debug = 0;                                                     # debug variable (0: normal, 1: debug mode)

# defining the URLs and commands for queries
my $arrtypus = 'arr';                                              # define 'arr' as the type for arrivals
my $deptypus = 'dep';                                              # define 'dep' as the type for departures
my $arrurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$arrtypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $depurl = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$deptypus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $arrcommand = "elinks -dump -dump-width 200 \"$arrurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data
my $depcommand = "elinks -dump -dump-width 200 \"$depurl\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

my $rheight=10;
my $rwidth=5;

# call the main subroutine
checkNet ($num);                                                   # call the main subroutine once at the beginning

############
# end main #
############

#
# subroutine to run the checking command, parse its output and 
# call the update subroutine
#

sub checkNet {
  my $num = shift;                                                        # define input variables
  my @text;                                                               # declare final arrays for shell (@text)                          
  my @deptext; my @arrtext;                                               # declare intermediate arrays (shell output) for departures (@deptext) and arrivals (@arrtext)
  open (ACPI, "$arrcommand |") || die "can't open pipe!";                 # open filestream ACPI reading the output from arrcommand
  Fahrinfo_ubuntu::read_finfo(\*ACPI,$station,\@arrtext);                 # subroutine to read data from ACPI and write it to @arrtext
  close ACPI;                                                             # close filestream ACPI
  open (ACPJ, "$depcommand |") || die "can't open pipe!";                 # open filestream ACPJ reading the output from depcommand
  Fahrinfo_ubuntu::read_finfo(\*ACPJ,$station,\@deptext);                 # subroutine to read data from ACPJ and write it to @deptext
  close ACPJ;                                                             # close filestream ACPJ
  my $natext = $#arrtext;                                                 # define size of @arrtext
  my $ndtext = $#deptext;                                                 # define size of @deptext
  Fahrinfo_ubuntu::pmax_dual(\@arrtext,\@deptext);                        # cut size of @arrtext, @deptext
  Fahrinfo_ubuntu::cmp_dual(\@arrtext,\@deptext,\@text);                  # compose data from @arrtext and @deptext to @text
  my $nrtext = $#text;                                                    # define size of @text
  if($nrtext > -1){                                                       # if @text is not empty
    print "@text";                                                        # output to the shell
  }
}

#
# subroutine to find suitable stations/bus stops to the input string and output them as a list (if more than one) or return correct station/bus stop name and number to the main program
#
sub haltnummer {
  my $halt = shift;                                                                               # read input string
  open(DATA, "$FindBin::RealBin/../data/fahrinfo-elinks2.dat") || die "can't open 'fahrinfo-elinks2.dat'";          # open filestream reading the local file fahrinfo-elinks2.dat
  my $numma = -1;                                                                                 # number of highest entry in @haltarr and @haltnrarr
  my $str;                                                                                        # UTF-8 decoded version of lines in fahrinfo-elinks2.dat
  my $haltu;                                                                                      # UTF-8 decoded version of $halt
  my @haltarr;                                                                                    # final array of station/bus stop names
  my @haltnrarr;                                                                                  # final array of station/bus stop numbers
  while (my $line = <DATA>) {                                                                     # read from fahrinfo-elinks2.dat
    $haltu = decode_utf8($halt);                                                                  # UTF-8 decode $halt to $haltu
    $str = decode_utf8($line);                                                                    # UTF-8 decode $line to $str
    if(index($line, $halt) != -1) {                                                               # if strings are identical
      $numma++;                                                                                   # increase number of highest entry
      $str =~ /(\"[a-zA-Z0-9\.\,\-\/\+\[\]\(\) ]+\")([\t])([0-9]+)/;                              # split line to name and number of station/bus stop
      push(@haltarr,$1);                                                                          # write name of station/bus stop to @haltarr
      push(@haltnrarr,$3);                                                                        # write number of station/bus stop to @haltnrarr
    }
  }
  my @testarr;                                                                                    # return array
  if($numma == 0){                                                                                # if only one station/bus stop is found
    @testarr = ($haltarr[$numma],$haltnrarr[$numma]);                                             # return array contains the name and the number of the chosen station/bus stop
    return @testarr;                                                                              # return the array to the main program
  }elsif($numma > 0 and defined $xnumm){                                                          # if more than one station/bus stop is found, but the chosen one is selected by $xnumm
    @testarr = ($haltarr[$xnumm],$haltnrarr[$xnumm]);                                             # return array contains the name and the number of the chosen station/bus stop
    return @testarr;                                                                              # return the array to the main program
  }else{                                                                                          # otherwise (no station/bus stop chosen or found)
    my $lnum = $numma + 1;                                                                        # number of found stations/bus stops
    print "$lnum Haltestellen gefunden: \n";                                                      # print number of found stations/bus stops
    foreach my $axnum (0..$#haltarr) {                                                            # print every suitable station/bus stop with a preceding number
      print "$axnum : $haltarr[$axnum]\n";
    }
    exit;                                                                                         # exit the program
  }
}

sub print_exit {                                                                                  # subroutine defining the error text (wrong input) and exiting the program
  print "Gebrauch: finfo-dual.pl \"Haltestelle\" <code> <Zeile (optional)>\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}

