#!/usr/bin/perl -w

use strict;
use warnings;
use utf8;
use Switch;
use FindBin;

my $browser = join('',$ARGV[0]);                                                         # Browser, in dem die Seiten aufgerufen werden sollen (w3m oder elinks)
my $modus = join('',$ARGV[1]);                                                           # Modus, in dem die Seiten aufgerufen werden (mobil oder stat)
my $listfile = join('',$ARGV[2]);                                                        # Datei, aus der die Adressen entnommen werden
my $mdate = `date | sed 's/ä/a/g'`;                                                      # speichere Ausgabe von date in $mdate
my @madate = split_date($mdate);                                                         # Aufspaltung in einzelne Bestandteile, Speicherung im @madate
my $fullthread = agent_thread(@madate);                                                  # Bestimmung, wie viel Daten zu bearbeiten sind (0: Mo-Fr außerhalb der Kernzeit, 1: Mo-Fr in der Kernzeit, 
                                                                                         # 2: wochenends/feiertags außerhalb der Kernzeit, 3: wochenends/feiertags in der Kernzeit)
my $number;
my @names;
my $nam;
open(LIST, $listfile) || die "Cannot open $listfile\n";
my $hol = 0;
if ($fullthread >= 4) {
 $hol = 1;
}else{
 $hol = 0;
}
my $innkernz = 0;
my $outkernz = 0;
if ($fullthread == 7 or $fullthread == 3){
 $innkernz = 1;
 $outkernz = 1;
}
if ($fullthread == 5 or $fullthread == 1){
 $outkernz = 1;
}
while ( my $line = <LIST> ){
   if ($line =~ /\#/){
     next;
   }
   if ($hol == 1){
     if ($line =~ /\%/){
       next;
     }
   }
   if ($innkernz == 0 and $outkernz == 1){
     if ($line =~ /\&/) {
       next;
     }
   }
   if ($innkernz == 0 and $outkernz == 0){
     if ($line =~ /\§/ or $line =~ /\&/){
       next;
     }
     if ($hol == 1 and $line =~ /\?/){
       next;
     }
   }
   if ( $line =~ /\%/ ){
    $nam = depercent($line);
   }elsif( $line =~ /\§/ ){
    $nam = deparagr($line);
   }elsif( $line =~ /\&/ ){
    $nam = deamp($line);
   }elsif( $line =~ /\?/ ){
    $nam = deques($line);
   }else{
    $nam = $line;
   }
   chomp($nam);
   push(@names,$nam);
}
my $shortmode = substr($modus,0,1);
system ("sh $FindBin::Bin/twitter.sh $browser $shortmode @names");

sub split_date {
 my $str = shift;
 $str =~ /([DFMS][aior]) ([0-9]{1,2})\. ([ADFJMNOS][aäekopu][bgilnprtvz]) ([0-9]{2})\:([0-9]{2})\:([0-9]{2}) ([\w]{3,4}) ([0-9]{4})/;
 my @arr;
 $arr[0] = $1;                                                                               # Wochentag
 $arr[1] = $2;                                                                               # Tag
 my $mon = $3;                                                                               # Monat (Wort)
 $arr[2] = wordtonumber_month($mon);                                                         # Monat (Zahl)
 $arr[3] = $4;                                                                               # Stunde
 $arr[4] = $5;                                                                               # Minute
 $arr[5] = $6;                                                                               # Sekunde
 $arr[6] = $7;                                                                               # Zeitzone
 $arr[7] = $8;                                                                               # Jahr
 return @arr;
}

sub wordtonumber_month {
     # input : an English or German abbreviation of a month, output: a number in the range [0:11]
    my $str = shift;
    $str =~/([a-zA-Z]{3})/;
    my $oldmonth = $1;
    my $newmonth;
    switch ($oldmonth) {
        case /Jan/      {$newmonth = 1;}
        case /Feb/      {$newmonth = 2;}
        case /Mar/      {$newmonth = 3;}
        case /Apr/      {$newmonth = 4;}
        case /Ma[iy]/   {$newmonth = 5;}
        case /Jun/      {$newmonth = 6;}
        case /Jul/      {$newmonth = 7;}
        case /Aug/      {$newmonth = 8;}
        case /Sep/      {$newmonth = 9;}
        case /O[ck]t/   {$newmonth = 10;}
        case /Nov/      {$newmonth = 11;}
        case /De[cz]/   {$newmonth = 12;}
    }
    return $newmonth;
}

sub calc_ifholiday {
    my @array = @_;
    my $weekday = $array[0];
    my $day = $array[1];
    my $month = $array[2];
    my $year = $array[7];
    my $holiday = 0;
    if ( $weekday =~ /Sa/ or $weekday =~ /So/ ){
       $holiday = 1;
       goto fin_calc_ifholiday;
    }else{
       $holiday = 0;
    }
    switch ($month){
        case [2,4,6,7,8,9,11] {$holiday = 0; goto fin_calc_ifholiday;}                # months without a holiday with fix date
        case [1,5]            {
                    switch ($day){
                        case [2..31] {$holiday = 0; goto fin_calc_ifholiday;}
                        case 1       {$holiday = 1; goto fin_calc_ifholiday;}
                    }
        }
        case 10               {
                    switch ($day){
                        case 3 {$holiday = 1; goto fin_calc_ifholiday;}
                        else   {$holiday = 0; goto fin_calc_ifholiday;}
                    }
        }
        case 12               {
                    switch ($day){
                        case [24..26] {$holiday = 1; goto fin_calc_ifholiday;}
                        case 31       {$holiday = 1; goto fin_calc_ifholiday;}
                        else          {$holiday = 0; goto fin_calc_ifholiday;}
                    }
        }
    }
    fin_calc_ifholiday:
    return $holiday;
}

sub calc_ifeaster {
    my @array = @_;
    my $weekday = $array[0];
    my $day = $array[1];
    my $month = $array[2];
    my $year = $array[7];
    my $holiday = 0;
    # calculus for Easter date
    my $k = ganz($year/100);                                                                  # secular number
    my $m = 15 + ganz((3*$k+3)/4) - ganz((8*$k+13)/25);                                       # secular moon equation
    my $s = 2 - ganz((3*$k+3)/4);                                                             # secular sun equation
    my $a = $year % 19;                                                                       # moon parameter
    my $d = (19*$a+$m) % 30;                                                                  # initial guess for first full moon in spring
    my $r = ganz(($d + ganz($a/11))/29);                                                      # correction for first full moon in spring
    my $og = 21 + $d - $r;                                                                    # Easter limit
    my $sz = 7 - ($year + ganz($year/4) + $s) % 7;                                            # first Sunday in March
    my $oe = 7 - ($og - $sz) % 7;                                                             # difference between first Sunday in March and Easter Sunday
    my $os = $og + $oe;                                                                       # Easter Sunday as a March date (also numbers larger than 31 are possible)
    my $kf = $os - 2;                                                                         # Good Friday as a March date
    my $om = $os + 1;                                                                         # Easter Monday as a March date
    my $ch = $os + 39;                                                                        # Ascension of Jesus Day as a March date
    my $pm = $os + 50;                                                                        # Pentecost Monday as a March date
    my @kfd;                                                                                  # "real" Good Friday date
    my @omd;                                                                                  # "real" Easter Monday date
    my @chd;                                                                                  # "real" Ascension of Jesus date
    my @pmd;                                                                                  # "real" Pentecost Monday date
    switch ($kf) {
        case [32..61] {$kfd[0] = $kf - 31; $kfd[1] = 4;}
        case [1..31]  {$kfd[0] = $kf;      $kfd[1] = 3;}
    }
    switch ($om) {
        case [32..61] {$omd[0] = $om - 31; $omd[1] = 4;}
        case [1..31]  {$omd[0] = $om;      $omd[1] = 3;}
    }
    switch ($ch) {
        case [93..122]{$chd[0] = $ch - 92; $chd[1] = 6;}
        case [62..92] {$chd[0] = $ch - 61; $chd[1] = 5;}
        case [32..61] {$chd[0] = $ch - 31; $chd[1] = 4;}
    }
    switch ($pm) {
        case [93..122]{$pmd[0] = $pm - 92; $pmd[1] = 6;}
        case [62..92] {$pmd[0] = $pm - 61; $pmd[1] = 5;}
        case [32..61] {$pmd[0] = $pm - 31; $pmd[1] = 4;}
    }
    if($day == $kfd[0] and $month == $kfd[1] and $weekday =~ /Fr/){
        $holiday = 1;
        goto fin_calc_ifeaster;
    }
    if($day == $omd[0] and $month == $omd[1] and $weekday =~ /Mo/){
        $holiday = 1;
        goto fin_calc_ifeaster;
    }
    if($day == $chd[0] and $month == $chd[1] and $weekday =~ /Do/){
        $holiday = 1;
        goto fin_calc_ifeaster;
    }
    if($day == $pmd[0] and $month == $pmd[1] and $weekday =~ /Mo/){
        $holiday = 1;
        goto fin_calc_ifeaster;
    }
    fin_calc_ifeaster:
    return $holiday;
}

sub ganz {
 my $str = shift;
 $str =~ /([0-9]+)[\.]*[0-9]*/;
 my $res = $1;
 return $res;
}

sub agent_thread {
 my @arr = @_;
 my $ifholiday = calc_ifholiday(@arr);
 my $ifeaster = calc_ifeaster(@arr);
 my $hora = $arr[3];
 my $min = $arr[4];
 my $res = 0;
 if ( $ifholiday == 1 or $ifeaster == 1 ){
   $res = 4;
 }
 if ($res == 4){
   if ($hora >= 9 and $hora <= 17){
     $res = 7;
     goto fin_agentthread;
   }elsif ($hora == 18 and $min == 0){
     $res = 7;
     goto fin_agentthread;
   }elsif ($hora >= 7 and $hora < 9){
     $res = 5;
     goto fin_agentthread;
   }elsif ($hora >= 18 and $hora <= 20){
     $res = 5;
     goto fin_agentthread;
   }elsif ($hora == 21 and $min == 0){
     $res = 5;
     goto fin_agentthread;
   }else{
     $res = 4;
     goto fin_agentthread;
   }
 }elsif($res == 0){
   if ($hora >= 7 and $hora <= 20){
     $res = 3;
     goto fin_agentthread;
   }elsif ($hora == 21 and $min == 0){
     $res = 3;
     goto fin_agentthread;
   }elsif ($hora == 6 or $hora == 21){
     $res = 1;
     goto fin_agentthread;
   }elsif ($hora == 22 and $min == 0){
     $res = 1;
     goto fin_agentthread;
   }else{
     $res = 0;
     goto fin_agentthread;
   }
 }
 fin_agentthread:
 return $res;
}

sub depercent {
 my $str = shift;
 $str =~ /\%([a-zA-Z\_\#0-9]*)/;
 my $res = $1;
 return $res;
}

sub deparagr {
 my $str = shift;
 $str =~ /\§([a-zA-Z\_\#0-9]*)/;
 my $res = $1;
 return $res;
}

sub deamp {
 my $str = shift;
 $str =~ /\&([a-zA-Z\_\#0-9]*)/;
 my $res = $1;
 return $res;
}

sub deques {
 my $str = shift;
 $str =~ /\?([a-zA-Z\_\#0-9]*)/;
 my $res = $1;
 return $res;
}
