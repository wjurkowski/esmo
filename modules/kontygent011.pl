#!/usr/bin/perl -w
##parametry programu:[lista plikow z sekwencjami] [katalog z sekwencjami] [szeroko?? okna]
#Wiktor Jurkowski. Please cite: Brylinski M.,Konieczny L.,Czerwonko P.,Jurkowski W. and Roterman I.
#Early-Stage Folding in Proteins (In Silico)Sequence-to-Structure Relation, J Biomed Biotechnol. 2005; 2005(2): 65–79.
#w kazdym plku sekwencja strukturalna i AA

use strict;
use warnings;
use defCT;
use defsubCT;
use CTcalc;

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if ($#ARGV != 0) {die "Program uzywany z parametrami!: [plik z parametrami]\n";}
open(OPCJE, "< $ARGV[0]") or die "Can not open an input file: $!";#plik z parametrami

#constans and predeclared variables
my @aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
my @sctab=('A','B','C','D','E','F','G');
my ($lin);

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
my $czas=localtime(time());
my $me =getlogin();	

my $wyniki="main.out";#wyniki kontrolne zapisane w pliku
$wyniki=$params{'out_file'} if ($params{'out_file'});
open(GOUT,">$wyniki") or die "Can not write an output file: m1 $!";

printf GOUT "Main ouput file generated with perl script: kontyngent.pl (W.Jurkowski) \n";
printf GOUT "User: %s Time: %s\n", $me, $czas;
printf GOUT "Selected parameters were used: \n";
while(my($key, $value)=each(%params)){
printf GOUT "$key = $value\n";}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#kazdy plik po kolei
if($params{'run_type'} eq 'define'){
my ($fname,@aasekw,@sasekw,$line);

my $inplist=$params{'inpseq_list'};
open(INPUT1, "< $inplist") or die "Can not open an input file: m2 $inplist $!";
my @lista=<INPUT1>;#wczytanie pliku wejsciowego, nazwy plikow z sekwencjami w tablicy @lista
close (INPUT1);
chomp @lista;

foreach $fname (@lista){#po plikach
print "chuj $fname\n";
my $sekw=$params{'inpseq_dir'};
#$katalog="$sekw";
#chdir "$katalog" or die "Nie moge wejsc do $sekw: $!\n";
open(INPUT2, "< $sekw/$fname") or die "Can’t open input file: m3 $sekw/$fname $!";
printf GOUT "file analyzed: $fname\n";

my @data=<INPUT2>;
close (INPUT2);
chomp @data;
my $nl=-1;
my $cc=0;
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
my $naa=$#aasekw;#dlugosc sekwencji AA
my $nsa=$#sasekw;#dlugosc sekwencji SA

if($params{'CTtype'} eq 'normseq'){#tworzenie AAvsSA
defseqCT(\%params,\@aasekw,\@sasekw,\$naa,\$fname);
}
elsif($params{'CTtype'} eq 'normstr'){#tworzenie SAvsAA
defstrCT(\%params,\@aasekw,\@sasekw,\$nsa,\$fname);
}

}#koniec plikow
}#koniec bloku dla definicji

elsif($params{'run_type'} eq 'analyze'){#analiza
printf GOUT "Programs run in analyze mode: Existing CT is beeing analyzed\n";
 if($params{'subbase'} eq 'limit'){#definicja podbaz
  if($params{'limit_row'} eq 'SAletters' and $params{'limit_col'} eq 'hydrophob'){
  
  defsubCT(\%params);}
 } 
 elsif($params{'subbase'} eq 'condseq'){#creates conditional AAvsSA
 defcondCT(\%params); 
 }#koniec definicji podbaz
 
 if($params{'inform'}==1){#analiza informacji

calSE(\%params);

 }#koniec analizy informacji	
}#koniec bloku dla analizy
my $koniec=localtime(time());
printf GOUT "Run completed. Time: $koniec\n";




