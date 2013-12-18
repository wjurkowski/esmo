# CondCT - defined constants for identifying codesets
# $Id: CondCT.pm,v 2.7 2004/06/10 21:19:34 neilb Exp $ 
package defcondCT;

use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw();


sub defcct{

$cct="condCT1.dat";
open(SOUT,"> $cct") or die "Can’t write output file: $!";#open the cct1 table

for $i (0..$naa-$wind){# dla kazdego AA przesuwamy okno	
	
	$s1="";#motyw strukturalny odpowiadajacy sekwencyjnemu
	for $k($i..$i+$wind){
	$s1=$s1.$sasekw[$k];}
	
	if ( -f "$a1") {#jezeli istnieje 
	$czy=0;
	
	my @plik=<SOUT>;
	chomp @plik;
	close (SOUT);
		foreach $row (@plik){
		 if($row=~/$s1/){
		  $czy=1;	
		 @wart=split /\t/,$row;
		 $wart[1]++; 
		 $row=join "\t", $wart[0],$wart[1];}	
		}
		if($czy==1){#i jest w nim odpowiedni motyw strukturalny
		open(SOUT,"> $a1") or die "I can’t open and write output file: $!";
		for $i (0..$#plik){
		print SOUT "$plik[$i]\n";}
		}
		if($czy==0){#jezeli istnieje ale niema takiego motywu strukturalnego
		open(SOUT,">> $a1") or die "I can’t open and write output file: $!";
		print SOUT "$s1\t1\n";}
	}
 }



}

































sub cct{
foreach $kat(@aatab){#motywy sekwencyjne
 	chdir "$aalibry/$kat/" or die "Nie moge wejsc do katalogu $aalibry/$kat $!\n";
	foreach $AA2(@aatab){
  	 foreach $AA3(@aatab){
	  foreach $AA4(@aatab){
		$motyw=$kat.$AA2.$AA3.$AA4;
		if( -f "$motyw"){#jezeli istnieje otworz go
		$cmn++;#numer kolumny
		open(CTAAINP,"< $motyw") or die "Can’t open SA motifs file of AA motif: $!, $motyw\n";
		my @plik=<CTAAINP>;#wczytaj do tablicy
		chomp @plik;
		close (CTAAINP);
#		 $rmn=0;
#		 $motsum[$cmn]=0;
#		 $SEC2[$cmn]=0;
#		 $SEC7[$cmn]=0;
#		 for $i (0..$#plik){#strukturalny dla danego sekwencyjnego
#		 $rmn++;#liczba stanow dla danej kolumny
#		 @para=split(/\s+/,$plik[$i]);
#		 $samot[$i]=$para[0];
#		 $rcnt[$i]=$para[1];
#		 $motsum[$cmn]=$motsum[$cmn]+$para[1];}#nastepny strukturalny
#		 $maxpr=0;

		 	
		 
		 	

		 for $i (0..$#plik){#prawdopodobienstwa
		 $pr[$i]=$rcnt[$i]/$motsum[$cmn];
			if($maxpr < $pr[$i]){
			$maxpr=$pr[$i];
			$maxprsa=$samot[$i];}
		 $entr2=$pr[$i]*(log($pr[$i])/log(2));
		 $entr7=$pr[$i]*(log($pr[$i])/log(7));
		 $SEC2[$cmn]=$SEC2[$cmn]+$entr2;
		 $SEC7[$cmn]=$SEC7[$cmn]+$entr7;
		 }
		 $MPSA[$cmn]=$maxpr;
		 if($rmn!=0){
		 $SEC2[$cmn]=-$SEC2[$cmn];
		 $SEC7[$cmn]=-$SEC7[$cmn];
		 $SECMAX2[$cmn]=log($rmn)/log(2);
		 $SECMAX7[$cmn]=log($rmn)/log(7);
		 $dSEC2[$cmn]=$SECMAX2[$cmn]-$SEC2[$cmn];
		 $dSEC7[$cmn]=$SECMAX7[$cmn]-$SEC7[$cmn];
		 printf SECOUT "%d\t%s\t%s\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n",$cmn,$motyw,$maxprsa,$SEC2[$cmn],$SECMAX2[$cmn],$dSEC2[$cmn],$SEC7[$cmn],$SECMAX7[$cmn],$dSEC7[$cmn];}#nastepny istniejacy
		 elsif($rmn==0){
		 print "WARNING: empty motif $motyw \n";}
		 }
	  }#nastepny motyw			
	 }
	}
 	}#nastepny katalog z motywami sekwencyjnymi	

 
}