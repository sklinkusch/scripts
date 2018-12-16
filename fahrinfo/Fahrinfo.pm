package Fahrinfo;
use strict;
use FindBin;
use Switch;
use List::Util qw(reduce any all none notall first max maxstr min minstr product sum sum0 pairs unpairs pairkeys pairvalues pairgrep pairfirst pairmap shuffle uniq uniqnum uniqstr);

sub add_linebreaks {
  my $in  = $_[0];
  my $ine = $_[1];
  foreach my $x (0..$#$in){
    $$in[$x] = sprintf("%s\n",$$in[$x]);
    $$ine[$x] = sprintf("%s\n",$$ine[$x]) if $x != 0;
    $$ine[$x] = sprintf(" %s\n",$$ine[$x]) if $x == 0;
  }
}

sub calculate_differenz {
   my ($plantime, $realtime) = (@_);
   $plantime =~ /([012][0-9])\:([012345][0-9])/;
   my $planstd = $1;
   my $planmte = $2;
   $realtime =~ /([012][0-9])\:([012345][0-9])/;
   my $realstd = $1;
   my $realmte = $2;
   my $differenz;
   $differenz = 60*($realstd - $planstd) + ($realmte - $planmte);
   if ($differenz < -120){
     $differenz += 1440;
   }
   return $differenz;
}

sub calc_filter {
  my $filtra = shift;
  my $fbahn = 0; my $rbahn = 0; my $sbahn = 0; my $ubahn = 0; my $tram = 0; my $bus = 0; my $ferry = 0;
  RETURN:
  switch ($filtra) {
    case [64..127]    {$fbahn = 1; $filtra -= 64; goto RETURN; }
    case [32..63]     {$rbahn = 1; $filtra -= 32; goto RETURN; }
    case [16..31]     {$sbahn = 1; $filtra -= 16; goto RETURN; }
    case [8..15]      {$ubahn = 1; $filtra -= 8;  goto RETURN; }
    case [4..7]       {$tram = 1; $filtra -= 4;   goto RETURN; }
    case [2..3]       {$bus = 1; $filtra -= 2;    goto RETURN; }
    case [1]          {$ferry = 1;}
  }
  my $filtnum = sprintf("%1d%1d%1d%1d%1d%1d%1d%1d%1d", $sbahn, $ubahn, $tram, $bus, $ferry, $fbahn, $rbahn, 0, 1);
  return $filtnum;
}

sub calc_hour {
  my $ctime = shift;
  $ctime = /([\d]{2})\:([\d]{2})/;
  my $hour = $1;
  return $hour;
}

sub compose_dual {
  my $ina = $_[0];
  my @inax = @{$_[0]};
  my $inea = $_[1];
  my @ineax = @{$_[1]};
  my $ind = $_[2];
  my @indx = @{$_[2]};
  my $ined = $_[3];
  my @inedx = @{$_[3]};
  my $t = $_[4];
  my $et = $_[5];
  my $nat = $#inax;
  my $ndt = $#indx;
  my $nz = min($nat,$ndt);
  foreach my $xa (0..$nz){
    if($xa == 0){
      $$t[$xa] = $inax[$xa];
      $$et[$xa] = " " . $ineax[$xa];
#      $$et[$xa] = $ineax[$xa];
    }else{
      $$t[$xa] = $inax[$xa];
      $$et[$xa] = $ineax[$xa];
    }
  }
  if($ndt > $nat){
    foreach my $xz (($nz+1)..$ndt){
      $$t[$xz] = sprintf("%-85s",$main::spaced);
      $$et[$xz] = sprintf("%-85s",$main::spaced);
    }
  }elsif($ndt < $nat){
    foreach my $xz (($nz+1)..$nat){
      push(@indx,"\n");
      push(@inedx,"\n");
    }
  }
  foreach my $xd (0..$ndt){
    $$t[$xd] = $$t[$xd] . $indx[$xd] . "\n";
    $$et[$xd] = $$et[$xd] . $inedx[$xd] . "\n";
  }
}

sub compose_fstrict_dual {
  my $ina = $_[0];
  my @inax = @{$_[0]};
  my $inea = $_[1];
  my @ineax = @{$_[1]};
  my $ind = $_[2];
  my @indx = @{$_[2]};
  my $ined = $_[3];
  my @inedx = @{$_[3]};
  my $t = $_[4];
  my $et = $_[5];
  my $nat = $#inax;
  my $ndt = $#indx;
  my $spaced = ' ';
  foreach my $xa (0..$nat){
    if($xa == 0){
      $$t[$xa] = sprintf("%-87s",$inax[$xa]);
      $$et[$xa] = sprintf(" %-87s",$ineax[$xa]);
#      $$et[$xa] = $ineax[$xa];
    }else{
      $$t[$xa] = $inax[$xa];
      $$et[$xa] = $ineax[$xa];
    }
  }
  if($ndt > $nat){
    foreach my $xz (($nat+1)..$ndt){
      $$t[$xz] = sprintf("%-87s",$spaced);
      $$et[$xz] = sprintf("%-87s",$spaced);
    }
  }elsif($ndt < $nat){
    foreach my $xz (($ndt+1)..$nat){
      push(@indx,"\n");
      push(@inedx,"\n");
    }
  }
  foreach my $xd (0..$ndt){
    $$t[$xd] = $$t[$xd] . $indx[$xd] . "\n";
    $$et[$xd] = $$et[$xd] . $inedx[$xd] . "\n";
  }
}

sub def_spec {
  my ($dtime,$A,$d,$D,$f,$F,$s,$S,$N) = @_;
  my $spec;
  switch ($dtime){
    case {$dtime <= $A}          {$spec = 'A';}
    case [($A+1)..-1]            {$spec = 'a';}
    case [0]                     {$spec = 'P';}
    case [1..$d]                 {$spec = 'd';}
    case [($d+1)..$D]            {$spec = 'D';}
    case [($D+1)..$f]            {$spec = 'f';}
    case [($f+1)..$F]            {$spec = 'F';}
    case [($F+1)..$s]            {$spec = 's';}
    case [($s+1)..$S]            {$spec = 'S';}
    case [($S+1)..$N]            {$spec = 'N';}
    case {$dtime > $N}           {$spec = 'x';}
    else                         {$spec = '-';}
  }
  return $spec;
}

sub define_product {
   my @arr = @_;
   my $linie = $arr[0];
   my $fahrtziel = $arr[1];
   my $spec = '0';
   $spec = 'F' if ($linie =~ /ICE|IC|EC|RJ|D|EN|NJ|TGV|THA|HBX|FLX/);
   $spec = 'R' if ($linie =~ /[I]?[R]{1}[BE]{1}/);
   $spec = 'S' if ($linie =~ /[S]{1}(1|2|25|26|3|41|42|45|46|47|5|7|75|8|85|9)$/);
   $spec = 'S' if ($linie =~ /^[0-9]{4,6}$/ and $fahrtziel =~ /^S /);
   $spec = 'S' if ($linie =~ /^[0-9]{4,6}$/ and $fahrtziel =~ /^Ringbahn S4[12]/);
   $spec = 'U' if ($linie =~ /[U]{1}(1|2|3|4|5|55|6|7|8|9)$/);
   $spec = 'M' if ($linie =~ /^[M]{1}(1|2|4|5|6|8|10|13|17)$/);
   $spec = 'm' if ($linie =~ /^[M]{1}(11|19|21|27|32|37|41|44|45|46|48|49|76|77|82|85)$/);
   $spec = 'T' if ($linie =~ /^[1-9]{1}[0-9]{1}$/);
   $spec = 'x' if ($linie =~ /^[X]{1}[1-9][0-9]?$/);
   $spec = 'x' if ($linie =~ /TXL/);
   $spec = 'b' if ($linie =~ /^[1-9]{1}[0-9]{2}$/);
   $spec = 'n' if ($linie =~ /^[N]{1}[1-9][0-9]?$/);
   $spec = 'f' if ($linie =~ /^[F]{1}[0-9]{1,2}$/);
   $spec = 'e' if ($spec eq '0');
   return $spec;
}

sub define_specifier {
   my @inp = @_;
   my $dtime = $inp[0];
   my $tspec = $inp[1];
   my $ctime = $inp[2];
   my $hour = Fahrinfo::calc_hour($ctime);
   my $spec;
   switch ($tspec) {
     case 'b'            {goto BUS;}
     case 'e'            {goto ELSE;}
     case 'F'            {goto FERN;}
     case 'f'            {goto FERRY;}
     case 'M'            {switch ($hour) {
                            case [0..4]     {goto NACHT;}
                            case [5..23]    {goto METRO;}
                            else            {goto METRO;}
                         }}
     case 'm'            {switch ($hour) {
                            case [0..4]     {goto NACHT;}
                            case [5..23]    {goto METRO;}
                            else            {goto METRO;}
                         }}
     case 'n'            {goto NACHT;}
     case 'R'            {goto REGIO;}
     case 'S'            {goto SBAHN;}
     case 'T'            {goto TRAM;}
     case 'U'            {goto UBAHN;}
     case 'x'            {goto BUS;}
   }
   FERN:
   $spec = Fahrinfo::def_spec($dtime,-10,9,29,59,89,119,179,209);
   goto AUSW;
   REGIO:
   $spec = Fahrinfo::def_spec($dtime,-5,4,9,19,29,59,89,119);
   goto AUSW;
   SBAHN:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,4,9,14,19,39);
   goto AUSW;
   UBAHN:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,3,5,7,9,19);
   goto AUSW;
   METRO:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,3,5,7,9,29);
   goto AUSW;
   TRAM:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,4,9,14,19,39);
   goto AUSW;
   BUS:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,4,9,14,19,39);
   goto AUSW;
   NACHT:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,4,9,19,29,59);
   goto AUSW;
   FERRY:
   $spec = Fahrinfo::def_spec($dtime,-2,1,2,3,4,5,6,9);
   goto AUSW;
   ELSE:
   $spec = Fahrinfo::def_spec($dtime,-3,2,4,7,9,14,19,39);
   goto AUSW;
   AUSW:
   return $spec;
}

sub filter_sgl {
  my $fls = $_[0];
  my $fl  = $_[1];
  my $ot  = $_[2];
  my $oet = $_[3];
  my $ft  = $_[4];
  my $fet = $_[5];
  my $and = 0;
  my @open;
  my @close;
  my @ort = @$ot;
  my @oret = @$oet;
  my @prt;
  my @pret;
  push (@$ft,$$ot[0]);
  push (@$fet,$$oet[0]);
  $and = 1 if (join(' ',@$fl) =~ /[\{\}]/);
  if ($and == 1){
    foreach my $xz (0..$#$fl){
      if ($$fl[$xz] =~ /^\{.+$/){
        push(@open,$xz);
        $$fl[$xz] = substr($$fl[$xz],1);
      }
      if ($$fl[$xz] =~ /^.+\}$/){
        push(@close,$xz);
        $/ = "\}";
        chomp($$fl[$xz]);
      }
    }
    main::print_exit() if(($#open - $#close) != 0);
    foreach my $xz (0..$#open){
      main::print_exit() if(($close[$xz] - $open[$xz]) != 1);
    }
    if ($fls eq 'f'){
      foreach my $xz (1..$#$ot){
        foreach my $xy (0..$#open){
          if ($$ot[$xz] =~ /$$fl[$open[$xy]]/ and $$ot[$xz] =~ /$$fl[$close[$xy]]/){
            push(@prt,$$ot[$xz]);
            push(@pret,$$oet[$xz]);
            last;
          }
        }
      }
    }else{
      foreach my $xz (1..$#$ot){
        foreach my $xy (0..$#open){
          if ($$ot[$xz] =~ /$$fl[$open[$xy]]/ and $$ot[$xz] =~ /$$fl[$close[$xy]]/){
            $$ot[$xz]  = undef;
            $$oet[$xz] = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (0..$#open){
      $$fl[$open[$xy]]   = undef;
      $$fl[$close[$xy]] = undef;
    }
  }
  if ($fls eq 'f'){
    foreach my $xz (1..$#$ot) {
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($$ot[$xz] =~ /$$fl[$xy]/){
            push(@prt,$$ot[$xz]);
            push(@pret,$$oet[$xz]);
            last;
          }
        }
      }
    }
  }else{
    foreach my $xz (1..$#$ot) {
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($$ot[$xz] =~ /$$fl[$xy]/){
            $$ot[$xz]  = undef;
            $$oet[$xy] = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (1..$#$ot) {
      if (defined $$ot[$xy]) {
        push(@prt,$$ot[$xy]);
        push(@pret,$$oet[$xy]);
      }
    }
  }
  foreach my $xy (1..$#ort){
    foreach my $xz (0..$#prt){
      if ($ort[$xy] eq $prt[$xz]){
        push(@$ft,$ort[$xy]);
        push(@$fet,$oret[$xy]);
        last;
      }
    }
  }
}

sub filter_sgl_p {
  my $fls = $_[0];
  my $fl  = $_[1];
  my $ot  = $_[2];
  my $oet = $_[3];
  my $ft  = $_[4];
  my $fet = $_[5];
  my $and = 0;
  my $dumstring;
  my @open;
  my @close;
  my @ort = @$ot;
  my @oret = @$oet;
  my @prt;
  my @pret;
  $and = 1 if (join(' ',@$fl) =~ /[\{\}]/);
  if ($and == 1){
    foreach my $xz (0..$#$fl){
      if ($$fl[$xz] =~ /^\{.+$/){
        push(@open,$xz);
        $$fl[$xz] = substr($$fl[$xz],1);
      }
      if ($$fl[$xz] =~ /^.+\}$/){
        push(@close,$xz);
        $/ = "\}";
        chomp($$fl[$xz]);
      }
    }
    main::print_exit() if(($#open - $#close) != 0);
    foreach my $xz (0..$#open){
      main::print_exit() if(($close[$xz] - $open[$xz]) != 1);
    }
    if ($fls eq 'f'){
      foreach my $xz (0..$#$ot){
        foreach my $xy (0..$#open){
          $dumstring = substr($$ot[$xz],21,64);
          if ($dumstring =~ /$$fl[$open[$xy]]/ and $dumstring =~ /$$fl[$close[$xy]]/){
            push(@prt,$$ot[$xz]);
            push(@pret,$$oet[$xz]);
            last;
          }
        }
      }
    }else{
      foreach my $xz (0..$#$ot){
        foreach my $xy (0..$#open){
          $dumstring = substr($$ot[$xz],21,64);
          if ($dumstring =~ /$$fl[$open[$xy]]/ and $dumstring =~ /$$fl[$close[$xy]]/){
            $$ot[$xz]  = undef;
            $$oet[$xz] = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (0..$#open){
      $$fl[$open[$xy]]   = undef;
      $$fl[$close[$xy]] = undef;
    }
  }
  if ($fls eq 'f'){
    foreach my $xz (0..$#$ot) {
      $dumstring = substr($$ot[$xz],21,64);
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($dumstring =~ /$$fl[$xy]/){
            push(@prt,$$ot[$xz]);
            push(@pret,$$oet[$xz]);
            last;
          }
        }
      }
    }
  }else{
    foreach my $xz (0..$#$ot) {
      $dumstring = substr($$ot[$xz],21,64);
      foreach my $xy (0..$#$fl) {
        if (defined $$fl[$xy]){
          if ($dumstring =~ /$$fl[$xy]/){
            $$ot[$xz]  = undef;
            $$oet[$xy] = undef;
            last;
          }
        }
      }
    }
    foreach my $xy (0..$#$ot) {
      if (defined $$ot[$xy]) {
        push(@prt,$$ot[$xy]);
        push(@pret,$$oet[$xy]);
      }
    }
  }
  foreach my $xy (0..$#ort){
    foreach my $xz (0..$#prt){
      if ($ort[$xy] eq $prt[$xz]){
        push(@$ft,$ort[$xy]);
        push(@$fet,$oret[$xy]);
        last;
      }
    }
  }
}

sub get_number {
 my $str = shift;
 $str =~ /([0-9]{1})([0-9]{6})/;
 my $add = '00';
 my $res = $1 . $add . $2;
 return $res;
}

sub make_fstrict {
  my $pt  = $_[0];
  my $pet = $_[1];
  my $t   = $_[2];
  my $et  = $_[3];
  my $exs;
  my $ini;
  my $ine;
  foreach my $x (0..$#$pt){
    if($x != 0){
      my $exs = substr($$pt[$x],87,55);
      if($exs =~ /^[ ]*$/ or $exs =~ /Gleis/){
        $ini = substr($$pt[$x],0,87);
        push(@$t,$ini);
        $ine = substr($$pet[$x],0,87);
        push(@$et,$ine);
      }else{
        next;
      }
    }else{
      $ini = substr($$pt[$x],0,87);
      push(@$t,$ini);
      $ine = substr($$pt[$x],0,87);
      push(@$et,$ine);
    }
  }
}

sub popmax_dual {
  my $ina = $_[0];
  my $inea = $_[1];
  my $ind = $_[2];
  my $ined = $_[3];
  my $nat = $#$ina;
  my $ndt = $#$ind;
  my $max = 35;
  my $diff;
  my $tempvar;
  if($nat > $max){
    $diff = $nat - $max;
    while ($diff > 0){
      $tempvar = pop(@$ina);
      $tempvar = pop(@$inea);
      $diff--;
    }
  }
  if($ndt > $max){
    $diff = $ndt - $max;
    while ($diff > 0){
      $tempvar = pop(@$ind);
      $tempvar = pop(@$ined);
      $diff--;
    }
  }
}

sub popmax_sgl {
  my $in  = $_[0];
  my $ine = $_[1];
  my $nt  = $#$in;
  my $max = 35;
  my $diff;
  my $tempvar;
  if($nt > $max){
    $diff = $nt - $max;
    while ($diff > 0){
      $tempvar = pop(@$in);
      $tempvar = pop(@$ine);
      $diff--;
    }
  }
}

sub read_fahrinfo {
  my $fs = $_[0];
  my $station = $_[1];
  my $t = $_[2];
  my $et = $_[3];
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $product;
  my $specifier;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $encline;
  my $spaced = ' ';
  my $fahrt;
  while (my $line = <$fs>){
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf("%-86s",$begline);
      $encline = $begaline;
      push(@$t,$begaline);
      push(@$et,$encline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            if(defined $realzeit){
              $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
              $specifier = Fahrinfo::define_specifier($difftime,$product,$currtime);
            }else{
              $difftime = 0;
              $specifier = '?';
            }
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z0-9]{2,3})[ ]{1,}/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([ ]{1})([A-Z0-9]{2,3})([ ]{1,})/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $skip_request == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = 'X';
            my $resttime = Fahrinfo::calculate_differenz($planzeit, $currtime);
            $difftime = undef;
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $skip_request == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and defined $planzeit and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ \t]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            if(defined $realzeit){
              $difftime = Fahrinfo::calculate_differenz($planzeit,$realzeit);
            }else{
              $difftime = 0;
            }
            $specifier = 'Z';
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            if (defined $difftime){
              my $fziel = Fahrinfo::remove_blanks($fahrtziel);
              if (defined $realzeit){
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel);
              }else{
                if (defined $linie){
                  $fahrt = sprintf("%-5s ----- %-1s %+4d $-1s %-9s %-55s",$planzeit,$specifier,$difftime,$product,$linie,$fziel);
                }else{
                  $fahrt = sprintf("%-5s ----- %-1s %+4d %-1s %-9s %-55s",$planzeit,$specifier,$difftime,$product,$spaced,$fziel);
                }
              }
            }else{
              my $fziel = Fahrinfo::remove_blanks($fahrtziel);
              if (defined $realzeit){
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel);
              }else{
                if (defined $linie){
                  $fahrt = sprintf("%-5s ----- %-1s %-4s %-1s %-9s %-55s",$planzeit,$specifier,$spaced,$product,$linie,$fziel);
                }else{
                  $fahrt = sprintf("%-5s ----- %-1s %-4s %-1s %-9s %-55s",$planzeit,$specifier,$spaced,$product,$spaced,$fziel);
                }
              }
            }
            push(@$t,$fahrt);
            push(@$et,$fahrt);
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
}

sub read_params {
  my $delt = $_[0];
  open (PARAM, "$FindBin::Bin/Fahrinfo.param") || die "cannot open params file\n";
  while (my $line = <PARAM>){
    if($line =~ /^[0-9]{1,}/){
      $line =~ /^([0-9]{1,})/;
      $delt = $1;
      last;
    }
  }
  close PARAM;
  return $delt;
}

sub read_strict {
  my $fs = $_[0];
  my $station = $_[1];
  my $t = $_[2];
  my $et = $_[3];
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $def_extra = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $specifier;
  my $extra_stop;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $encline;
  my $product;
  my $spaced = ' ';
  while (my $line = <$fs>) {
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf ("%-86s",$begline);
      $encline = $begaline;
      push(@$t,$begline);
      push(@$et,$encline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9][ ]{1,}[A-Z]{0,3}[ ]?[0-9]{0,5}[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis [A-Za-z0-9\,\.\-\+\/\(\)\[\]]*/){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $extra_stop = $6;
            $def_extra = 1;
            $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $1;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $2;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $difftime = undef;
            $specifier = '?';
            $extra_stop = $3;
            $def_extra = 1;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $specifier = 'X';
            $difftime = undef;
            my $resttime = Fahrinfo::calculate_differenz($planzeit, $currtime);
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            $specifier = 'Z';
            $difftime = Fahrinfo::calculate_differenz($planzeit,$realzeit);
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $curline =~ /Haltestelle\/Gleis/) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
            $curline =~ /Haltestelle\/Gleis ([A-Za-z0-9\+\-\.\,\/\(\)\[\] ]*)/;
            $extra_stop = $1;
            $def_extra = 1;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            my $fahrt;
            if ($def_extra == 1 and defined $extra_stop){
              if (defined $difftime){
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }
            }else{
              if (defined $difftime){
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$spaced) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$spaced) unless ($specifier eq 'x');
              }
            }
            push(@$t,$fahrt) unless ($specifier eq 'x');
            push(@$et,$fahrt) unless ($specifier eq 'x');
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            $def_extra = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
}

sub read_strict_p {
  my $fs = $_[0];
  my $station = $_[1];
  my $t = $_[2];
  my $et = $_[3];
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $def_extra = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $specifier;
  my $extra_stop;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $encline;
  my $product;
  my $spaced = ' ';
  while (my $line = <$fs>) {
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf ("%-86s",$begline);
      $encline = $begaline;
#      push($t,$begline);
#      push($et,$encline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9][ ]{1,}[A-Z]{0,3}[ ]?[0-9]{0,5}[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis [A-Za-z0-9\,\.\-\+\/\(\)\[\]]*/){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $extra_stop = $6;
            $def_extra = 1;
            $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $1;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $2;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $difftime = undef;
            $specifier = '?';
            $extra_stop = $3;
            $def_extra = 1;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $specifier = 'X';
            $difftime = undef;
            my $resttime = Fahrinfo::calculate_differenz($planzeit, $currtime);
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            $specifier = 'Z';
            $difftime = Fahrinfo::calculate_differenz($planzeit,$realzeit);
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $curline =~ /Haltestelle\/Gleis/) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
            $curline =~ /Haltestelle\/Gleis ([A-Za-z0-9\+\-\.\,\/\(\)\[\] ]*)/;
            $extra_stop = $1;
            $def_extra = 1;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            my $fahrt;
            if ($def_extra == 1 and defined $extra_stop){
              if (defined $difftime){
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }
            }else{
              if (defined $difftime){
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$station) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$station) unless ($specifier eq 'x');
              }
            }
            if ($fahrt =~ /$station/){
              push(@$t,$fahrt) unless ($specifier eq 'x');
              push(@$et,$fahrt) unless ($specifier eq 'x');
            }
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            $def_extra = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
}

sub read_strict_pn {
  my $fs = $_[0];
  my $first = $_[1];
  my $station = $_[2];
  my $t = $_[3];
  my $et = $_[4];
  my @pt;
  my @pet;
  my $index = 0;
  my $def_fahrtziel = 0;
  my $def_planzeit = 0;
  my $def_realzeit = 0;
  my $def_linie = 0;
  my $def_extra = 0;
  my $skip_request = 0;
  my $fahrtziel;
  my $planzeit;
  my $realzeit;
  my $linie;
  my $specifier;
  my $extra_stop;
  my $difftime;
  my $currtime;
  my $begline;
  my $begaline;
  my $encline;
  my $product;
  my $spaced = ' ';
  while (my $line = <$fs>) {
    if ( $line !~ /Abfahrt \(/ and $line !~ /Ankunft \(/ and $index == 0) {
      next;
    }
    if ( $line =~ /Abfahrt \(/ or $line =~ /Ankunft \(/){
      $line =~ /(A[bn][fk][au][hn][rf]t) \(([012][0-9]\:[012345][0-9])/;
      my $word = $1;
      $currtime = $2;
      $begline = sprintf ("%-.55s %-8s %-5s",$station,$word,$currtime);
      $begaline = sprintf ("%-86s",$begline);
      $encline = $begaline;
#      push($t,$begline);
#      push($et,$encline);
      $index = 1;
      next;
    }
    if ( $index == 1 and $line !~ /References/) {
      my $curline = Fahrinfo::subst_text($line);
      if ($curline =~ /In Richtung/ or $curline =~ /Aus Richtung/) {
        next;
      }
      if ($curline =~ /[ ]{10,24}[A-Za-z\,\-\.0-9\/\+\(\)\[\] ]{8,65}/ and $curline !~ /Zwischenhalte/ and $curline !~ /[012][0-9]\:[012345][0-9]/ and $curline !~ / (Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez) / and $curline !~ /taeglich/ and $curline !~ / (Mo|Di|Mi|Do|Fr|Sa|So) / and $curline !~ /Haltestelle\/Gleis/){ # and $def_fahrtziel == 0 and $def_linie == 0 and $def_planzeit == 0 and $def_realzeit == 0) {
        $curline =~ /([ ]{10,24})([A-Za-z\,\.\-0-9\/\+\(\)\[\] ]*)/;
        $fahrtziel = $2;
        $def_fahrtziel = 1;
        next;
      }
      if ($def_fahrtziel == 1){
        if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9]/){
          if ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9][ ]{1,}[A-Z]{0,3}[ ]?[0-9]{0,5}[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis [A-Za-z0-9\,\.\-\+\/\(\)\[\]]*/){
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $extra_stop = $6;
            $def_extra = 1;
            $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}[012][0-9]\:[012345][0-9] [012][0-9]\:[012345][0-9]/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{1,5})([012][0-9]\:[012345][0-9])[ ]{1}([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = $3;
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $difftime = Fahrinfo::calculate_differenz($planzeit, $realzeit);
            $specifier = Fahrinfo::define_specifier($difftime, $product,$currtime);
          }elsif ($curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /^[ ]{1,5}([012][0-9]\:[012345][0-9])[ ]{1}([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}[012][0-9]\:[012345][0-9][ ]*[A-Za-z0-9\,\.\-\/\+\(\)\[\] ]*Haltestelle\/Gleis ([A-Za-z0-9\,\.\-\+\/\(\)\[\] ]*)/;
            $planzeit = $1;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $2;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $difftime = undef;
            $specifier = '?';
            $extra_stop = $3;
            $def_extra = 1;
          }elsif ($curline =~ /Ausfall/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])[ ]{1}(Ausfall)([ ]*)([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]*[012][0-9]\:[012345][0-9]/;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = 'XXXXX';
            $def_realzeit = 1;
            $linie = $5;
            $def_linie = 1;
            $specifier = 'X';
            $difftime = undef;
            my $resttime = Fahrinfo::calculate_differenz($planzeit, $currtime);
#            if ($resttime >= 30){
#              $def_planzeit = 0;
#              $def_realzeit = 0;
#              $def_linie = 0;
#              $def_fahrtziel = 0;
#              next;
#            } 
          }elsif ($curline =~ /Zusatzfahrt/ and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1})(Zusatzfahrt)([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1,}([012][0-9]\:[012345][0-9])/;
            $planzeit = $2;
            $def_planzeit = 1;
            $skip_request = 1;
            if (defined $6 and length($6) > 0){
              $linie = $6;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            next;
          }elsif ($skip_request == 1 and $def_planzeit == 1 and $def_realzeit == 0){
            $curline =~ /([ ]{3})([012][0-9]\:[012345][0-9])([ ]{1,})([A-Z]{0,3}[ ]?[0-9]{0,5})[ ]{1}[012][0-9]\:[012345][0-9]/;
            $realzeit = $2;
            $def_realzeit = 1;
            if ($def_linie == 0){
              $linie = $4;
              $def_linie = 1;
              $product = Fahrinfo::define_product($linie,$fahrtziel);
            }
            $specifier = 'Z';
            $difftime = Fahrinfo::calculate_differenz($planzeit,$realzeit);
            $skip_request = 0;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0 and $curline =~ /Haltestelle\/Gleis/) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
            $curline =~ /Haltestelle\/Gleis ([A-Za-z0-9\+\-\.\,\/\(\)\[\] ]*)/;
            $extra_stop = $1;
            $def_extra = 1;
          }elsif ($skip_request == 0 and $def_planzeit == 0 and $def_realzeit == 0 and $def_linie == 0) {
            $curline =~ /^([ ]{1,5})([012][0-9]\:[012345][0-9])([\s]*)([A-Z]{0,3}[ ]?[0-9]{0,5}) /;
            $planzeit = $2;
            $def_planzeit = 1;
            $realzeit = '-----';
            $def_realzeit = 1;
            $linie = $4;
            $def_linie = 1;
            $product = Fahrinfo::define_product($linie,$fahrtziel);
            $specifier = '?';
            $difftime = undef;
          }
          if ($def_planzeit == 1 and $def_realzeit == 1 and $def_linie == 1) {
            my $fahrt;
            if ($def_extra == 1 and defined $extra_stop){
              if (defined $difftime){
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$extra_stop) unless ($specifier eq 'x');
              }
            }else{
              if (defined $difftime){
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %+4d %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$difftime,$product,$linie,$fziel,$station) unless ($specifier eq 'x');
              }else{
                my $fziel = Fahrinfo::remove_blanks($fahrtziel);
                $fahrt = sprintf("%-5s %-5s %-1s %-4s %-1s %-9s %-55s %-55s",$planzeit,$realzeit,$specifier,$spaced,$product,$linie,$fziel,$station) unless ($specifier eq 'x');
              }
            }
            if ($fahrt =~ /$station/){
              push(@pt,$fahrt) unless ($specifier eq 'x');
              push(@pet,$fahrt) unless ($specifier eq 'x');
            }
            $def_planzeit = 0;
            $def_realzeit = 0;
            $def_linie = 0;
            $def_fahrtziel = 0;
            $def_extra = 0;
            next;
          }
        }
      }
    }
    if ( $index == 1 and $line =~ /References/) {
      last;
    }
  }
  my $ref_pt = \@pt;
  my $ref_pet = \@pet;
  push(@$t,$ref_pt);
  push(@$et,$ref_pet);
}

sub remove_blanks {
 my $str = shift;
 $str =~ s/^\s*//;
 return $str;
}

sub search_linie {
  my $curline = shift;
  my $linie;
  if($curline =~ / [A-Z]{2,3} [0-9]{1,5} /){
    $curline =~ / ([A-Z]{2,3} [0-9]{1,5}) /;
    $linie = $1;
  }elsif($curline =~ /TXL/){
    $linie = 'TXL';
  }elsif($curline =~ / [MNSUX]{1}[0-9]{1,2} /){
    $curline =~ / ([MNSUX]{1}[0-9]{1,2}) /;
    $linie = $1;
  }elsif($curline =~ / [0-9]{3} /){
    $curline =~ / ([0-9]{3}) /;
    $linie = $1;
  }
  return $linie;
}

sub sort_entries {
  my $at  = @_[0];
  my $aet = @_[1];
  my $bt  = @_[2];
  my $bet = @_[3];
  my $t   = @_[4];
  my $et  = @_[5];
  for (my $x = 0; $x <= ($#$at + $#$bt); $x++){
    if(defined $$at[0] and defined $$bt[0]){
      if($$at[0] lt $$bt[0]){
        push(@$t,$$at[0]);
        push(@$et,$$aet[0]);
        shift(@$at);
        shift(@$aet);
      }else{
        push(@$t,$$bt[0]);
        push(@$et,$$bet[0]);
        shift(@$bt);
        shift(@$bet);
      }
    }elsif(defined $$at[0]){
      push(@$t,$$at[0]);
      push(@$et,$$aet[0]);
      shift(@$at);
      shift(@$aet);
    }else{
      push(@$t,$$bt[0]);
      push(@$et,$$bet[0]);
      shift(@$bt);
      shift(@$bet);
    }
  }
}

sub sort_entries_n {
  my $nrtot = @_[0];                                                                               # total number of stations
  my $mj    = @_[1];                                                                               # maximum connections per station
  my $pt    = @_[2];                                                                               # twodimensional array containing unsorted connections (terminal output)
  my $pet   = @_[3];                                                                               # twodimensional array containing unsorted connections (perl-Tk output)
  my $t     = @_[4];                                                                               # final onedimensional array for sorted connections (terminal output)
  my $et    = @_[5];                                                                               # final onedimensional array for sorted connections (perl-Tk output)
  my @testarr;                                                                                     # test array (terminal output)
  my @testarre;                                                                                    # test array (perl-Tk output)
  my $date = `date`;                                                                               # print output of date command in variable
  my @splitdate = split(' ',$date);                                                                # split date into parts
  my $timefull = $splitdate[3];                                                                    # get variable for time (XX:XX:XX)
  my @timecomp = split(':',$timefull);                                                             # split time into parts
  my $time = sprintf("%02u:%02u",$timecomp[0],$timecomp[1]);                                       # write time (XX:XX) into variable
  my @testarr_s;                                                                                   # sorted test array (terminal output)
  my @testarre_s;                                                                                  # sorted test array (perl-Tk output)
  my @ind;                                                                                         # index array
  # set index array values initially to 0
  foreach my $nh (0..$nrtot){
    $ind[$nh] = 0;
  }
  # sorting procedure
  foreach my $nt (0..($nrtot*$mj)) {
  # write time into test arrays
#    $testarr[0] = $time;
#    $testarre[0] = $time;
  # push the first not already used connection for each station into test array
    foreach my $nr (0..($nrtot-1)){
      if (defined $$pt[$nr][(0+$ind[$nr])]){
        push(@testarr,$$pt[$nr][(0+$ind[$nr])]);
        push(@testarre,$$pet[$nr][(0+$ind[$nr])]);
      }else{
        push(@testarr,undef);
        push(@testarre,undef);
      }
    }
  # sort arrays (ASCII procedure)
    @testarr_s = sort {$a cmp $b} @testarr;
    @testarre_s = sort {$a cmp $b} @testarre;
    if ($time =~ /^2/ and $testarr_s[0] =~ /^0/ and $testarr_s[$#testarr_s] !~ /^0/){
      my $zz = $#testarr;
      while ($zz >= 0) {
        if($testarr_s[$zz] =~ /^0/){
          foreach my $zy (1..($#testarr-$zz)){
            if (defined $testarr[($zz+$zy)]){
              push(@$t,$testarr_s[($zz+$zy)]);
              push(@$et,$testarre_s[($zz+$zy)]);
              last;
            }
          }
          last;
        }
        $zz--;
      }
    }else{
      foreach my $zz (0..$#testarr){
        if (defined $testarr_s[$zz]){
          push(@$t,$testarr_s[$zz]);
          push(@$et,$testarre_s[$zz]);
          last;
        }
      }
    }
  # update index array
    foreach my $nz (0..$#testarr){
      if($testarr[$nz] eq $$t[$#$t]){
        $ind[$nz] = $ind[$nz] + 1;
        last;
      }
    }
  # remove elements from test arrays
    foreach my $nz (0..$#testarr){
      shift(@testarr);
      shift(@testarre);
    }
  }
}

sub subst_text {
      my $lina = shift;
      $lina =~ s/\[[1-9][0-9]\]/    /g;
      $lina =~ s/\[[1-9][0-9]{2}\]/     /g;
      $lina =~ s/\[[1-9][0-9]{3}\]/      /g;
      $lina =~ s/\[IMG\]/     /g;
      $lina =~ s/^[ \t]*$//g;
      return $lina;
}

1;
