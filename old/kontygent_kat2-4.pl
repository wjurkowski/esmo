#!/usr/bin/perl -w
##parametry programu:[lista plikow z sekwencjami] [katalog z sekwencjami] [szerokość okna]
#w kazdym plku sekwencja strukturalna i AA
#
#
#PODSTAWOWE

if ($#ARGV != 2) {die "Program uzywany z parametrami!\
        [lista plikow z sekwencjami] [katalog z sekwencjami] [szerokość okna]\n";}
#katalog roboczy
$rootdir=`pwd`;
chomp $rootdir;
#$rootdir=$root."/";
#print "chuj $rootdir\n";

#otwiera plik wejsciowy
open(INPUT1, "< $ARGV[0]") or die "Can not open an input file: $!";

#wczytanie pliku wejsciowego, nazwy plikow z sekwencjami w tablicy @lista
my @lista=<INPUT1>;
close (INPUT1);
chomp @lista;

#wyniki kontrolne zapisane w pliku 
$wyniki=$ARGV[0].".out";
open(GOUT,"> $rootdir/$wyniki") or die "Can not write an output file: $!";

#PARAMETRY

#szerokosc okna do obliczen tablicy kontyngencji
$window=$ARGV[2]-1;

@aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
@sctab=('A','B','C','D','E','F','G');

#definicja sciezek, tworzenie katalogow
$library='library-0507';
$libseq='sek-str';
$libstr='str-sek';

$libaadir="$rootdir/$library/$libseq";
$libscdir="$rootdir/$library/$libstr";

mkdir("$rootdir/$library", 0755) if (! -d "$rootdir/$library");
mkdir("$libaadir", 0755) if (! -d "$libaadir");
mkdir("$libscdir", 0755) if (! -d "$libscdir");
foreach $aa(@aatab){
mkdir("$libaadir/$aa", 0755) if (! -d "$libaadir/$aa");
#chdir "$libaadir/$aa" or die "Nie moge wejsc do $aa: $!\n";
}
foreach $sc(@sctab){
mkdir("$libscdir/$sc", 0755) if (! -d "$libscdir/$sc");
#chdir "$rootdir/$library/$libseq/$sc" or die "Nie moge wejsc do $sc: $!\n";
}

#GLOWNY
#kazdy plik po kolei

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
for $i (0..$naa-$window){
chdir "$rootdir/$library/$libseq/$aasekw[$i]/"; 
#	$k=$i+$window-1;
	$a1="";
	for $k ($i..$i+$window){
#	print "$aasekw[$k]\n";
	$a1=$a1.$aasekw[$k];
	}
#	print "chuj $a1\n";
	
#motyw strukturalny odpowiadajacy sekwencyjnemu
	$s1="";
	for $k ($i..$i+$window){
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



