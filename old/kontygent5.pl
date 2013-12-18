#!/usr/bin/perl -w
##parametry programu:[lista plikow z sekwencjami] [katalog z sekwencjami] [szeroko?? okna]
#w kazdym plku sekwencja strukturalna i AA
#
#
#PODSTAWOWE

if ($#ARGV != 2) {die "Program uzywany z parametrami!\
        [lista plikow z sekwencjami] [katalog z sekwencjami] [plik z parametrami]\n";}
#katalog roboczy
$rootdir=`pwd`;
chomp $rootdir;
#$rootdir=$root."/";
#print "chuj $rootdir\n";

#otwiera pliki wejsciowe
open(INPUT1, "< $ARGV[0]") or die "Can not open an input file: $!";
open(OPCJE, "< $ARGV[2]") or die "Can not open an input file: $!";

#wczytanie pliku wejsciowego, nazwy plikow z sekwencjami w tablicy @lista
my @lista=<INPUT1>;
close (INPUT1);
chomp @lista;

#wyniki kontrolne zapisane w pliku 
$wyniki=$ARGV[0].".out";
open(GOUT,"> $rootdir/$wyniki") or die "Can not write an output file: $!";

#PARAMETRY

#wczytanie opcji 
my @param=<OPCJE>;
close (OPCJE);
chomp @param;
my %params = ();
foreach $lin(@param){
@para=split(/\s+/,$lin);
$params{$para[0]}=$para[1];
}
print GOUT "parameters read correctly\n";

#szerokosc okna do obliczen tablicy kontyngencji
$wind=$params{'window'};

@aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
@sctab=('A','B','C','D','E','F','G');
@hydrofob=('A','V','L','I','M','F','W','Y','P','C');

$library=$params{'libdir'};

if($params{'run_type'} eq 'sequential' or $params{'run_type'} eq 'structural'){ 
mkdir("$library", 0755) if (! -d "$library");
}
if($params{'run_type'} eq 'sequential'){
 foreach $aa(@aatab){
 mkdir("$library/$aa", 0755) if (! -d "$library/$aa");
 #chdir "$libaadir/$aa" or die "Nie moge wejsc do $aa: $!\n";
 }	
}
elsif($params{'run_type'} eq 'structural'){
 foreach $sc(@sctab){
 mkdir("$library/$sc", 0755) if (! -d "$library/$sc");
 #chdir "$rootdir/$library/$libseq/$sc" or die "Nie moge wejsc do $sc: $!\n";
 }
} 

#GLOWNY
#kazdy plik po kolei
#if($params{'run_type'} ne 'analyze'){
if($params{'run_type'} eq 'sequential' or $params{'run_type'} eq 'structural'){
foreach $fname (@lista) {
$katalog="$rootdir/$ARGV[1]";
#print "$katalog\n";
chdir "$katalog" or die "Nie moge wejsc do $ARGV[1]: $!\n";
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
	 goto READ
	 }
	 else{
	 next LINIA
	 }
READ:		if($cc==1){
		@aasekw=split //,$data[$nl+1];
		}
		if($cc==2){
		@sasekw=split //,$data[$nl+1];
		} 
	}

#dlugosc sekwencji
$naa=$#aasekw;
#@smcount=();
$czy=0;
# dla kazdego AA przesuwamy okno
for $i (0..$naa-$wind){
chdir "$rootdir/$library/$aasekw[$i]/"; 
#	$k=$i+$wind-1;
	$a1="";
	for $k($i..$i+$wind){
#	print "$aasekw[$k]\n";
	$a1=$a1.$aasekw[$k];
	}
#	print "chuj $a1\n";
	
#motyw strukturalny odpowiadajacy sekwencyjnemu
	$s1="";
	for $k($i..$i+$wind){
	$s1=$s1.$sasekw[$k];
	}
#	print "dupa $s1\n";
#	$s1=$sasekw[$i].$sasekw[$i+1].$sasekw[$i+2].$sasekw[$i+3];
#	$smcount[$hit][$key[1]][$key[2]][$key[3]][$key[4]]++;

#jezeli taki plik jeszcze nie istnieje
	if (! -f "$a1") {
	open(SOUT,"> $a1") or die "Can’t write output file: $!";
	printf SOUT "$s1\t1\n";
	close (SOUT);
	next;
	}
#jezeli istnieje i jest w nim odpowiedni motyw strukturalny
	if ( -f "$a1") {
	$czy=0;
	open(SOUT,"< $a1") or die "Can’t write output file: $!";
	my @plik=<SOUT>;
	chomp @plik;
	close (SOUT);
#	print "@plik\n";
		foreach $row (@plik){
		 if($row=~/$s1/){
		  $czy=1;	
		 @wart=split /\t/,$row;
		 $wart[1]++; 
#		print "$wart[0]\t$wart[1]\n";
		 $row=join "\t", $wart[0],$wart[1];
		print "$row\n";
		 }	
		}
#zapisanie do pliku	
		if($czy==1){
		open(SOUT,"> $a1") or die "I can’t open and write output file: $!";
		for $i (0..$#plik){
		print SOUT "$plik[$i]\n";
		}
		}

#jezeli istnieje ale niema takiego motywu strukturalnego
		if($czy==0){
#		print "dupa\n";
		open(SOUT,">> $a1") or die "I can’t open and write output file: $!";
		print SOUT "$s1\t1\n";
		}
	}
	
}
}
}

elsif($params{'run_type'} eq 'analyze'){

	if($params{'limit_row'} eq 'SAletters'){
		if($params{'limit_col'} eq 'hydrophob'){
		my @CT=();
		 foreach $pre(@sctab){
		 #tworz plik dla kazdej litery osobno
		 open(TOUT, ">$pre._SA_vs_hydropbocity-CT.dat") or die "I can’t open and write output file: $!\n";
		 }
		}
	}

 foreach $kat(@aatab){
 chdir "$rootdir/$library/$kat/" or die "Nie moge wejsc do katalogu $kat $!\n";

	foreach $AA2(@aatab){
  	 foreach $AA3(@aatab){
	  foreach $AA4(@aatab){
		$hfobc=0;
		$motyw=$kat.$AA2.$AA3.$AA4;
		#jezeli istnieje taki plik
		if( -f "$motyw"){
		 #sprawdzaj hydrofobowosc motywu sekwencyjnego
 		 foreach $hfba(@hydrofob){
		  if($kat eq $hfba){$hfobc++;}
		  if($AA2 eq $hfba){$hfobc++;}
		  if($AA3 eq $hfba){$hfobc++;}
		  if($AA4 eq $hfba){$hfobc++;}	
		 }
		#otworz go
		open(CTINP,"< $motyw") or die "Can’t open SA motifs file of AA motif: $!, $motyw\n";
		#wczytaj do tablicy
		my @plik=<CTINP>;
		chomp @plik;
		close (CTINP);

		 #pokroj wiersze
		 for $i (0..$#plik){
		 my %seven=();
		 @para=split(/\s+/,$plik[$i]);
			#przepusc dany motyw strukturalny

#			foreach $slett($sctab){
#			$zlicz=~/$slett/g
#			}

			@string=split("",$para[0]);
#		print "kaka, @string\n";
			for $k (0..$wind-1){
			 foreach $slett(@sctab){
#		print "kutas,$slett,$string[$k]\n";
			  if($string[$k] eq $slett){
			  $seven{$slett}++;
#			print "dupa,$seven{$slett}\n";
			  }
			 }
			}
			#sprawdz klucze i wartosci
			$j=-1;
			while(($klucz,$wart)=each %seven){
#			while($wart=values %seven){
			$j++;
#			print "chuj $klucz,$wart,$hfobc,$j,$para[1]\n";
			$CT[$wart][$hfobc][$j]=$CT[$wart][$hfobc][$j]+$para[1];
			}
		 #nastepny strukturalny dla danego sekwencyjnego
		 }
		#nastepny istniejacy	
		}
	  #nastepny motyw
	  }			
	 }
	}
#nastepny katalog z motywami sekwencyjnymi
 }
#zapisanie tabeli kontyngencji
for $j(0..6){
	for $w(0..4){
printf TOUT ("%5I\t%5I\t%5I\t%5I\t%5I\n",$CT[$w][0][$j],$CT[$w][1][$j],$CT[$w][2][$j],$CT[$w][3][$j],$CT[$w][4][$j]);
	}
}

}



