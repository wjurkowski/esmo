# CondCT - defined constants for identifying codesets
# $Id: CondCT.pm,v 2.7 2004/06/10 21:19:34 neilb Exp $ 
package defCondCT;

use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(@chftab);

sub cct{

printf GOUT "Programs run in analyze mode: Existing CT is beeing analyzed\n";
$aalibry="$rootdir/$params{'libaadir'}";
$salibry="$rootdir/$params{'libsadir'}";
   if($params{'limit_row'} eq 'SAletters'){
   if($params{'limit_col'} eq 'hydrophob'){
#definicja tablicy (referencje!!!)
		for $j (0..6){
		 for $wart (0..4){
		  for $hfobc (0..4){
		  $CT[$j][$wart][$hfobc]=0;}
		 }
		}	
	foreach $kat(@aatab){
 	chdir "$aalibry/$kat/" or die "Nie moge wejsc do katalogu $kat $!\n";
	foreach $AA2(@aatab){
  	 foreach $AA3(@aatab){
	  foreach $AA4(@aatab){
		$hfobc=0;
		$motyw=$kat.$AA2.$AA3.$AA4;
		if( -f "$motyw"){#jezeli istnieje taki plik
 		 foreach $hfba(@hydrofob){#sprawdzaj hydrofobowosc motywu sekwencyjnego
		  if($kat eq $hfba){$hfobc++;}
		  if($AA2 eq $hfba){$hfobc++;}
		  if($AA3 eq $hfba){$hfobc++;}
		  if($AA4 eq $hfba){$hfobc++;}	
		 }
		open(CTINP,"< $motyw") or die "Can’t open SA motifs file of AA motif: $!, $motyw\n";#otworz go
		my @plik=<CTINP>;#wczytaj do tablicy
		chomp @plik;
		close (CTINP);
		 for $i (0..$#plik){#pokroj wiersze
		 my %seven=();
		 @para=split(/\s+/,$plik[$i]);
			@string=split("",$para[0]);#przepusc dany motyw strukturalny
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

	chdir "$rootdir" or die "Nie moge wejsc do katalogu $kat $!\n";#zapisanie tabeli kontyngencji
	for $j(0..6){#tworz plik dla kazdej litery osobno
	open(TOUT, ">$sctab[$j].SA_vs_hydroph-CT.dat") or die "I can’t open and write output file: $!\n";
	 for $w (0..4){
	 printf TOUT  ("%7d\t%7d\t%7d\t%7d\t%7d\n",$CT[$j][$w][0],$CT[$j][$w][1],$CT[$j][$w][2],$CT[$j][$w][3],$CT[$j][$w][4]);
	 #printf TOUT $CT[$j][$w][0],$CT[$j][$w][1],$CT[$j][$w][2],$CT[$j][$w][3],$CT[$j][$w][4];
	 }
	close (TOUT);
	}
   }	
  }
 
}