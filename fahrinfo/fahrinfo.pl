#!/usr/bin/perl 

#
# A Perl/Tk wrapper to graphically show the arrival and departure table for a certain station/tramway stop/bus stop
#
# It calls 'elinks', takes the content of the page, and writes it to # the graphical window. It is further updated every 10 seconds.
#
# Written by SK
#

use Tk;                                # use graphical window
use strict;                            # all variables have to be declared before use
use warnings;                          # print all errors and warnings to the shell
use utf8;                              # use unicode table for font encryption
use Encode;                            # packages to reencode text
use open ':encoding(utf8)';                                                      
use open ':std';                                                                 
use FindBin;                           # package to find an extra file for modules
use lib $FindBin::RealBin;             # absolute path where extra file is found
use Fahrinfo;                          # name of extra file (without .pm ending)

### set variables
if ($#ARGV != 2 and $#ARGV != 3) {
  print_exit();
}

my $haltestelle = join('',$ARGV[0]);
my $xnumm = join('',$ARGV[3]) if $#ARGV == 3;
my @haltestellenarr = haltnummer($haltestelle);
my $haltestellennum = $haltestellenarr[1];
my $haltestellenr = Fahrinfo::get_number($haltestellennum);
my $station = $haltestellenarr[0];

my $typus = join('',$ARGV[1]);
if ($typus ne 'dep' and $typus ne 'arr') {
  print_exit();
}

my $filtre = join('',$ARGV[2]);
$filtre = 127 if ($filtre < 1 or $filtre > 127);
my $filter = Fahrinfo::calc_filter($filtre);

my $num;
my $debug = 0;

my $del;                                                # delay time (seconds)
$del = Fahrinfo::read_params(\$del);                    # read parameter from external file
my $DELAY = $del * 1000; 				# delay time (milliseconds)

my $url = "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestellenr&boardType=$typus&maxJourneys=200&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&pageViewMode=PRINT&dirINPUT=&";
my $command = "elinks -dump -dump-width 200 \"$url\" | sed \'s/ä/ae/g;s/ö/oe/g;s/ü/ue/g;s/Ä/Ae/g;s/Ö/Oe/g;s/Ü/Ue/g;s/ß/ss/g;s/é/e/g\'";		# command to get fahrinfo data

my $rheight=10;
my $rwidth=5;

# create main window
my $mw = MainWindow->new();
$mw->geometry('800x550');
$mw->resizable(1, 0);

# set up quit key bindings
$mw->bind('<q>', sub {exit});
$mw->bind('<Escape>', sub {exit;});


my $label= $mw->Label(-font=>"Courier -12 bold",-justify=>'left')->pack;


# call the power checker now
checkNet ($label, $num);
$mw->repeat($DELAY, \&checkNet, $label, $num);

MainLoop;

############
# end main #
############

#
# subroutine to run the checking command, parse its output and 
# call the update subroutine
#
sub checkNet {
  my ($label, $num) = @_;
  my @text;  my @enctext;
  open (ACPI, "$command |") || die "can't open pipe!";
  Fahrinfo::read_fahrinfo(\*ACPI,$station,\@text,\@enctext);
  close ACPI;
  Fahrinfo::popmax_sgl(\@text,\@enctext);
  Fahrinfo::add_linebreaks(\@text,\@enctext);
  my $nrtext = $#text;
  if($nrtext > -1){
    print "@text" if ($debug);
    $label->configure(-text=>"@enctext");
  }
}

sub haltnummer {
  my $halt = shift;
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
  }elsif($numma > 0 and defined $xnumm){
    print "$xnumm : $haltarr[$xnumm]\n";
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
  print "Gebrauch: fahrinfo.pl \"Haltestelle\" <dep/arr> <code> <Zeile (optional)>\n";
  print "dep: Abfahrtsplan, arr: Ankunftsplan\n";
  print "code: Summe von Zahlen: Fernbahn (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)\n";
  exit;
}
