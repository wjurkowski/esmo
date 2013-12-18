#!/usr/bin/perl -w
##parametry programu:[lista plikow z sekwencjami] [katalog z sekwencjami] [szeroko?? okna]
#Wiktor Jurkowski. Please cite: Brylinski M.,Konieczny L.,Czerwonko P.,Jurkowski W. and Roterman I.
#Early-Stage Folding in Proteins (In Silico)Sequence-to-Structure Relation, J Biomed Biotechnol. 2005; 2005(2): 65–79.
#w kazdym plku sekwencja strukturalna i AA

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if ($#ARGV != 0) {die "Program uzywany z parametrami!: [plik z parametrami]\n";}
$rootdir=`pwd`;#katalog roboczy
chomp $rootdir;
open(OPCJE, "< $ARGV[0]") or die "Can not open an input file: $!";#plik z parametrami

#PARAMETRY
#wczytanie opcji 
my @param=<OPCJE>;
close (OPCJE);
chomp @param;
my %params = ();
my @para=();
foreach $lin(@param){
@para=split(/\s+/,$lin);
$params{$para[0]}=$para[1];
}
$czas=localtime(time());
$me =getlogin();	

$inplist=$params{'inpseq_list'};
open(INPUT1, "< $inplist") or die "Can not open an input file: $!";
my @lista=<INPUT1>;#wczytanie pliku wejsciowego, nazwy plikow z sekwencjami w tablicy @lista
close (INPUT1);
chomp @lista;
$wyniki=$inplist.".out";#wyniki kontrolne zapisane w pliku
open(GOUT,"> $rootdir/$wyniki") or die "Can not write an output file: $!";

printf GOUT "Main ouput file generated with perl script: kontyngent.pl (W.Jurkowski) \n";
printf GOUT "User: %s Time: %s\n", $me, $czas;
printf GOUT "Selected parameters were used: \n";
while(my($key, $value)=each(%params)){
printf GOUT "$key = $value\n";}

$wind=$params{'window'}-1;#szerokosc okna do obliczen tablicy kontyngencji
@aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
@sctab=('A','B','C','D','E','F','G');
@hydrofob=('A','V','L','I','M','F','W','Y','P','C');

if($params{'run_type'} eq 'sequential'){
 $aalibry="$rootdir/$params{'libaadir'}";
 mkdir("$aalibry", 0755) if (! -d "$aalibry");
 chdir "$aalibry" or die "Nie moge wejsc do $aalibry: $!\n";
 foreach $aa(@aatab){
 mkdir("$aalibry/$aa", 0755) if (! -d "$aalibry/$aa");}	
}
elsif($params{'run_type'} eq 'structural'){
 $salibry="$rootdir/$params{'libsadir'}";
 mkdir("$salibry", 0755) if (! -d "$salibry");
 chdir "$salibry" or die "Nie moge wejsc do $salibry: $!\n";
 foreach $sc(@sctab){
 mkdir("$salibry/$sc", 0755) if (! -d "$salibry/$sc");}
} 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#kazdy plik po kolei
if($params{'run_type'} eq 'sequential' or $params{'run_type'} eq 'structural'){
foreach $fname (@lista){#po plikach
$sekw=$params{'inpseq_dir'};
$katalog="$rootdir/$sekw";
chdir "$katalog" or die "Nie moge wejsc do $sekw: $!\n";
open(INPUT2, "< $fname") or die "Can’t open input file: $!";
printf GOUT "file analyzed: $fname\n";

my @data=<INPUT2>;
close (INPUT2);
chomp @data;
$nl=-1;
$cc=0;
LINIA:	foreach $line (@data){
	$nl++;
	 if ($line=~/^>/){
	 $cc++;
	 goto READ}
	 else{
	 next LINIA}
READ:		if($cc==1){
		@aasekw=split //,$data[$nl+1];}
		if($cc==2){
		@sasekw=split //,$data[$nl+1];} 
	}
$naa=$#aasekw;#dlugosc sekwencji AA
$nsa=$#sasekw;#dlugosc sekwencji SA
$czy=0;
if($params{'run_type'} eq 'sequential'){#tworzenie AAvsSA
 for $i (0..$naa-$wind){# dla kazdego AA przesuwamy okno	
	chdir "$aalibry/$aasekw[$i]/" or die "Nie moge wejsc do $aalibry/$aasekw[$i]/ $!\n";
	$a1="";#motyw sekwencyjny
	for $k($i..$i+$wind){
	$a1=$a1.$aasekw[$k];}
	$s1="";#motyw strukturalny odpowiadajacy sekwencyjnemu
	for $k($i..$i+$wind){
	$s1=$s1.$sasekw[$k];}
	if (! -f "$a1") {#jezeli taki plik jeszcze nie istnieje
	open(SOUT,"> $a1") or die "Can’t write output file: $!";
	printf SOUT "$s1\t1\n";
	close (SOUT);
	next;}
	if ( -f "$a1") {#jezeli istnieje 
	$czy=0;
	open(SOUT,"< $a1") or die "Can’t write output file: $!";
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
if($params{'run_type'} eq 'structural'){#tworzenie SAvsAA	
	for $i (0..$nsa-$wind){# dla kazdego AA przesuwamy okno	
	chdir "$salibry/$sasekw[$i]/" or die "Nie moge wejsc do $salibry/$aasekw[$i]/ $!\n";
	$a1="";#motyw sekwencyjny
	for $k($i..$i+$wind){
	$a1=$a1.$aasekw[$k];}
	$s1="";#motyw strukturalny odpowiadajacy sekwencyjnemu
	for $k($i..$i+$wind){
	$s1=$s1.$sasekw[$k];}
	if (! -f "$s1") {#jezeli taki plik jeszcze nie istnieje
	open(AOUT,"> $s1") or die "Can’t write output file: $s1 $fname $!";
	printf AOUT "$a1\t1\n";
	close (AOUT);
	next;}
	if ( -f "$s1") {#jezeli istnieje 
	$czy=0;
	open(AOUT,"< $s1") or die "Can’t write output file: $!";
	my @plik=<AOUT>;
	chomp @plik;
	close (AOUT);
		foreach $row (@plik){
		 if($row=~/$a1/){
		  $czy=1;	
		 @wart=split /\t/,$row;
		 $wart[1]++; 
		 $row=join "\t", $wart[0],$wart[1];}	
		}
		if($czy==1){#i jest w nim odpowiedni motyw strukturalny
		open(AOUT,"> $s1") or die "I can’t open and write output file: $!";
		for $i (0..$#plik){
		print AOUT "$plik[$i]\n";}
		}
		if($czy==0){#jezeli istnieje ale niema takiego motywu strukturalnego
		open(AOUT,">> $s1") or die "I can’t open and write output file: $!";
		print AOUT "$a1\t1\n";}
	}
 }
}
}#koniec plikow
}#koniec bloku dla definicji
elsif($params{'run_type'} eq 'analyze'){#analiza
printf GOUT "Programs run in analyze mode: Existing CT is beeing analyzed\n";
$aalibry="$rootdir/$params{'libaadir'}";
$salibry="$rootdir/$params{'libsadir'}";
 if($params{'subbase'}==1){#definicja podbaz
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
 }#koniec definicji podbaz
 if($params{'inform'}==1){#analiza informacji
$centrop=$params{'se_cout'};
$rentrop=$params{'se_rout'};
open(SECOUT,">$centrop") or die "Can’t write entropies for AA motifs (columns)\n";
open(SEROUT,">$rentrop") or die "Can’t write entropies for SA motifs (rows)\n";	
print SEROUT "row #\tmotif\tmax p\tSE2\tSEmax2\tdSE2\tSE20\tSEmax20\tdSE20\n";
print SECOUT "col #\tmotif\tmax p\tSE2\tSEmax2\tdSE2\tSE7\tSEmax7\tdSE7\n";
	$cmn=0;
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
		 $rmn=0;
		 $motsum[$cmn]=0;
		 $SEC2[$cmn]=0;
		 $SEC7[$cmn]=0;
		 for $i (0..$#plik){#strukturalny dla danego sekwencyjnego
		 $rmn++;#liczba stanow dla danej kolumny
		 @para=split(/\s+/,$plik[$i]);
		 $samot[$i]=$para[0];
		 $rcnt[$i]=$para[1];
		 $motsum[$cmn]=$motsum[$cmn]+$para[1];}#nastepny strukturalny
		 $maxpr=0;
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
	$rmn=0;
 foreach $kat(@sctab){#motywy strukturalne
 	chdir "$salibry/$kat/" or die "Nie moge wejsc do katalogu $kat $!\n";
	foreach $SA2(@sctab){
  	 foreach $SA3(@sctab){
	  foreach $SA4(@sctab){
		$motyw=$kat.$SA2.$SA3.$SA4;
		if( -f "$motyw"){#jezeli istnieje otworz go
		$rmn++;#numer wiersza
		open(CTSAINP,"< $motyw") or die "Can’t open SA motifs file of AA motif: $!, $motyw\n";
		my @plik2=<CTSAINP>;#wczytaj do tablicy
		chomp @plik2;
		close (CTSAINP);
		 $cmn=0;
		 $motsum[$rmn]=0;
		 $SER2[$rmn]=0;
		 $SER20[$rmn]=0;
		 for $i (0..$#plik2){#sekwencyjny dla danego strukturalnego
		 $cmn++;#liczba stanow dla danego wiersza
		 @para=split(/\s+/,$plik2[$i]);
		 $aamot[$i]=$para[0];
		 $ccnt[$i]=$para[1];
		 $motsum[$rmn]=$motsum[$rmn]+$para[1];}#nastepny sekwencyjny
		 $maxpr=0;
		 for $i (0..$#plik2){#prawdopodobienstwa
		 $pr[$i]=$ccnt[$i]/$motsum[$rmn];
			if($maxpr < $pr[$i]){
			$maxpr=$pr[$i];
			$maxpraa=$aamot[$i];}
		 $entr2=$pr[$i]*(log($pr[$i])/log(2));
		 $entr20=$pr[$i]*(log($pr[$i])/log(20));
		 $SER2[$rmn]=$SER2[$rmn]+$entr2;
		 $SER20[$rmn]=$SER20[$rmn]+$entr20;}
		 $MPAA[$rmn]=$maxpr;
		 if($cmn!=0){
		 $SER2[$rmn]=-$SER2[$rmn];
		 $SER20[$rmn]=-$SER20[$rmn];
		 $SERMAX2[$rmn]=log($cmn)/log(2);
		 $SERMAX20[$rmn]=log($cmn)/log(20);
		 $dSER2[$rmn]=$SERMAX2[$rmn]-$SER2[$rmn];
		 $dSER20[$rmn]=$SERMAX20[$rmn]-$SER20[$rmn];
		 printf SEROUT "%d\t%s\t%s\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n",$rmn,$motyw,$maxpraa,$SER2[$rmn],$SERMAX2[$rmn],$dSER2[$rmn],$SER20[$rmn],$SERMAX20[$rmn],$dSER20[$rmn];}#nastepny istniejacy
		 elsif($cmn==0){
		 print "WARNING: empty motif $motyw \n";}
		 }
	  }#nastepny motyw			
	 }
	}
 	}#nastepny katalog z motywami sekwencyjnymi
 }#koniec analizy informacji	
}#koniec bloku dla analizy
$koniec=localtime(time());
printf GOUT "Run completed. Time: $koniec\n";




