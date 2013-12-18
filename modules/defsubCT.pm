# defsubCT - defines subCT
# $Id: defsubCT.pm
package defsubCT;

use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(hfbvsSA);

sub hfbvsSA{
my @aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
my @sctab=('A','B','C','D','E','F','G');
my @hydrofob=('A','V','L','I','M','F','W','Y','P','C');
my($params)=@_;
my $wind=$$params{'window'}-1;#szerokosc okna do obliczen tablicy kontyngencji
my $rootdir=`pwd`;#katalog roboczy
my $aalibry="$$params{'libaadir'}";
my $salibry="$$params{'libsadir'}";
my (@CT,$kat,$AA2,$AA3,$AA4,$i,$j,$k,$w,$klucz,$wart,$hfobc,@seven,$hfba,$slett);

for $j (0..6){
 for $wart (0..4){
  for $hfobc (0..4){
  $CT[$j][$wart][$hfobc]=0;}
 }
}

chdir "$aalibry" or die "Nie moge wejsc do katalogu d1 $kat $!\n";
foreach $kat(@aatab){
 	chdir "$kat" or die "Nie moge wejsc do katalogu d2 $kat $!\n";
	foreach $AA2(@aatab){
  	 foreach $AA3(@aatab){
	  foreach $AA4(@aatab){
		$hfobc=0;
		my $motyw=$kat.$AA2.$AA3.$AA4;
		if( -f "$motyw"){#jezeli istnieje taki plik
 		 foreach $hfba(@hydrofob){#sprawdzaj hydrofobowosc motywu sekwencyjnego
		  if($kat eq $hfba){$hfobc++;}
		  if($AA2 eq $hfba){$hfobc++;}
		  if($AA3 eq $hfba){$hfobc++;}
		  if($AA4 eq $hfba){$hfobc++;}	
		 }
		open(CTINP,"< $motyw") or die "Can�t open SA motifs file of AA motif: $!, $motyw\n";#otworz go
		my @plik=<CTINP>;#wczytaj do tablicy
		chomp @plik;
		close (CTINP);
		 for $i (0..$#plik){#pokroj wiersze
		 my %seven=();
		 my @para=split(/\s+/,$plik[$i]);
			my @string=split("",$para[0]);#przepusc dany motyw strukturalny
			for $k (0..$wind-1){
			 foreach $slett(@sctab){
			  if($string[$k] eq $slett){
			  $seven{$slett}++;}
			 }
			}
			while(($klucz,$wart)=each %seven){#sprawdz klucze i wartosci
				for $j (0..6){
				 if($klucz eq $sctab[$j]){
				 $CT[$j][$wart][$hfobc]=$CT[$j][$wart][$hfobc]+$para[1];}
				}
			}
		 }#nastepny strukturalny dla danego sekwencyjnego
		}#nastepny istniejacy
	  }#nastepny motyw			
	 }
	}
 }#nastepny katalog z motywami sekwencyjnymi
chdir "$rootdir" or die "Nie moge wejsc do katalogu d3 $rootdir $!\n";#zapisanie tabeli kontyngencji
for $j(0..6){#tworz plik dla kazdej litery osobno
open(TOUT, ">$sctab[$j].SA_vs_hydroph-CT.dat") or die "I can�t open and write output file: $!\n";
 for $w (0..4){
 printf TOUT  ("%7d\t%7d\t%7d\t%7d\t%7d\n",$CT[$j][$w][0],$CT[$j][$w][1],$CT[$j][$w][2],$CT[$j][$w][3],$CT[$j][$w][4]);
 #printf TOUT $CT[$j][$w][0],$CT[$j][$w][1],$CT[$j][$w][2],$CT[$j][$w][3],$CT[$j][$w][4];
 }
close (TOUT);
} 
 
}#end of subroutine