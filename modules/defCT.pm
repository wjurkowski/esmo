# defCondCT - defines conditional CT
# $Id: defCondCT.pm
package defCT;

use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(defseqCT defstrCT defcondCT);

sub defseqCT{
my @aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
my @sctab=('A','B','C','D','E','F','G');
my($params,$aasekw,$sasekw,$naa,$fname)=@_;
my $wind=$$params{'window'}-1;#szerokosc okna do obliczen tablicy kontyngencji
my $czy=0;
my $rootdir=`pwd`;#katalog roboczy
chomp $rootdir;
my ($aa,$i,$k,$row,@wart);

my $aalibry="$$params{'libaadir'}";
 mkdir("$aalibry", 0755) if (! -d "$aalibry");
 chdir "$aalibry" or die "Nie moge wejsc do d1 $aalibry: $!\n";
 foreach $aa(@aatab){
 mkdir("$aa", 0755) if (! -d "$aa");}
 chdir "$rootdir" or die "Nie moge wejsc do d2: $rootdir $!\n";	

for $i (0..$$naa-$wind){# dla kazdego AA przesuwamy okno	
	chdir "$rootdir/$aalibry/$$aasekw[$i]/" or die "Nie moge wejsc do d3:  $rootdir/$aalibry/$$aasekw[$i]/ $!\n";
	my $a1="";#motyw sekwencyjny
	for $k($i..$i+$wind){
	$a1=$a1.$$aasekw[$k];}
	my $s1="";#motyw strukturalny odpowiadajacy sekwencyjnemu
	for $k($i..$i+$wind){
	$s1=$s1.$$sasekw[$k];}
	if (! -f "$a1") {#jezeli taki plik jeszcze nie istnieje
	open(SOUT,"> $a1") or die "Can?t write output file: $! $a1 for $fname";
	printf SOUT "$s1\t1\n";
	close (SOUT);
	next;}
	if ( -f "$a1") {#jezeli istnieje 
	$czy=0;
	open(SOUT,"< $a1") or die "Can?t write output file: $!";
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
		open(SOUT,"> $a1") or die "I can?t open and write output file: $!";
		for $i (0..$#plik){
		print SOUT "$plik[$i]\n";}
		}
		if($czy==0){#jezeli istnieje ale niema takiego motywu strukturalnego
		open(SOUT,">> $a1") or die "I can?t open and write output file: $!";
		print SOUT "$s1\t1\n";}
	}
 }
 chdir "$rootdir" or die "Nie moge wejsc do d4: $rootdir $!\n";	
 }#end of subroutine

sub defstrCT{
my @aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
my @sctab=('A','B','C','D','E','F','G');
my($params,$aasekw,$sasekw,$nsa,$fname)=@_;
my $wind=$$params{'window'}-1;#szerokosc okna do obliczen tablicy kontyngencji
my $czy=0;
my $rootdir=`pwd`;#katalog roboczy
chomp $rootdir;
my ($sc,$i,$k,$row,@wart);

my $salibry="$$params{'libsadir'}";
 mkdir("$salibry", 0755) if (! -d "$salibry");
 chdir "$salibry" or die "Nie moge wejsc do d5 $salibry: $!\n";
 foreach $sc(@sctab){
 mkdir("$sc", 0755) if (! -d "$sc");}
 chdir "rootdir" or die "Nie moge wejsc do d6 $rootdir: $!\n";

  	for $i (0..$$nsa-$wind){# dla kazdego AA przesuwamy okno	
 	chdir "$rootdir/$salibry/$$sasekw[$i]/" or die "Nie moge wejsc do d7 $salibry/$$aasekw[$i]/ $!\n";
 	my $a1="";#motyw sekwencyjny
 	for $k($i..$i+$wind){
 	$a1=$a1.$$aasekw[$k];}
 	my $s1="";#motyw strukturalny odpowiadajacy sekwencyjnemu
 	for $k($i..$i+$wind){
 	$s1=$s1.$$sasekw[$k];}
 	if (! -f "$s1") {#jezeli taki plik jeszcze nie istnieje
 	open(AOUT,"> $s1") or die "Can?t write output file: $s1 for $fname $!";
 	printf AOUT "$a1\t1\n";
 	close (AOUT);
 	next;}
 	if ( -f "$s1") {#jezeli istnieje 
 	$czy=0;
 	open(AOUT,"< $s1") or die "Can?t write output file: $!";
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
 		open(AOUT,"> $s1") or die "I can?t open and write output file: $!";
 		for $i (0..$#plik){
 		print AOUT "$plik[$i]\n";}
 		}
 		if($czy==0){#jezeli istnieje ale niema takiego motywu strukturalnego
 		open(AOUT,">> $s1") or die "I can?t open and write output file: $!";
 		print AOUT "$a1\t1\n";}
 	}
 }
chdir "$rootdir" or die "Nie moge wejsc do d8: $rootdir $!\n";	
}#end of subroutine

sub defcondCT{
my @aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
my @sctab=('A','B','C','D','E','F','G');
my $rootdir=`pwd`;#katalog roboczy
chomp $rootdir;
my($params)=@_;
my $aalibry="$$params{'libaadir'}";
my $caalibry="$$params{'cctaadir'}";
my ($aa,$kat,$AA2,$AA3,$AA4,$i,$j,$k,$row,@wart,@kod,@ilosc,$p4,$p3,$p2,$p1,$sa1,$sa2,$sa3);

 mkdir("$caalibry", 0755) if (! -d "$caalibry");
 chdir "$caalibry" or die "Nie moge wejsc do d9 $caalibry: $!\n";
 foreach $aa(@aatab){
 mkdir("$aa", 0755) if (! -d "$aa");}
 chdir "$rootdir" or die "Nie moge wejsc do d10 $rootdir: $!\n";

#chdir "$aalibry" or die "Nie moge wejsc do katalogu d11 $rootdir/$aalibry/$kat $!\n";
foreach $kat(@aatab){
# 	chdir "$rootdir/$aalibry/$kat/" or die "Nie moge wejsc do katalogu d11 $rootdir/$aalibry/$kat $!\n";
	foreach $AA2(@aatab){
  	 foreach $AA3(@aatab){
	  foreach $AA4(@aatab){
		my $motyw=$kat.$AA2.$AA3.$AA4;
		if( -f "$aalibry/$kat/$motyw"){#jezeli istnieje taki plik
 		open(CTINP,"< $aalibry/$kat/$motyw") or die "Can?t open SA motifs file of AA motif: $!, $aalibry/$kat/$motyw\n";#otworz go
		my @plik=<CTINP>;#wczytaj do tablicy
		chomp @plik;
		close (CTINP);
		open(CTOUT,"> $rootdir/$caalibry/$kat/$motyw") or die "Can?t open CCT file of AA motif: $!, $rootdir/$caalibry/$kat/$motyw\n";#print CCT
		
		 for $i (0..$#plik){#pokroj wiersze
		 my @para=split(/\s+/,$plik[$i]);
		 $kod[$i]=$para[0];
		 $ilosc[$i]=$para[1];
		 }
		 
		 my $sum0=0;
		 for $i (0..$#plik){#czworki
		 my $sum3=0;
		 $sum0=$sum0+$ilosc[$i];
		 my @str=split("",$kod[$i]);
		 my $mot3=$str[0].$str[1].$str[2];
		  for $j (0..$#plik){
		   if($kod[$j]=~/^$mot3/){
		   $sum3=$sum3+$ilosc[$j];}
		  }
		 $p4=$ilosc[$i]/$sum3;		  
		 printf CTOUT "%s\t%7.3f\n",$kod[$i],$p4;#print
		 }#nastepny strukturalny dla danego sekwencyjnego
  		  
 		 foreach $sa1(@sctab){#trojki
		  foreach $sa2(@sctab){
		  my $sum2=0;	
		  my $mot2=$sa1.$sa2;
		    for $j (0..$#plik){
		     if($kod[$j]=~/^$mot2/){
		     $sum2=$sum2+$ilosc[$j];}
		    }
		     foreach $sa3(@sctab){
		     my $sum3=0;	
    	   	     my $mot3=$sa1.$sa2.$sa3;
   		      for $j (0..$#plik){
		       if($kod[$j]=~/^$mot3/){
		       $sum3=$sum3+$ilosc[$j];}
		      }
			if($sum3 == 0){
			$p3=0;}
			elsif($sum3 > 0){
			$p3=$sum3/$sum2;
			printf CTOUT "%s\t%7.3f\n",$mot3,$p3;}#print
		     }
		   }
		 }
		 
		 foreach $sa1(@sctab){#dwojki
		 my $sum1=0;
		 my $mot1=$sa1;
		  for $j (0..$#plik){
		   if($kod[$j]=~/^$mot1/){
		   $sum1=$sum1+$ilosc[$j];}
		  }
		   foreach $sa2(@sctab){
		   my $sum2=0;	
		   my $mot2=$sa1.$sa2;
		    for $j (0..$#plik){
		     if($kod[$j]=~/^$mot2/){
		     $sum2=$sum2+$ilosc[$j];}
		    }
			if($sum2 == 0){
			$p2=0;}
			elsif($sum2>0){
			$p2=$sum2/$sum1;
			printf CTOUT "%s\t%7.3f\n",$mot2,$p2;}#print
		   }
		 } 
		 
		 foreach $sa1(@sctab){#jedynki
		 my $sum1=0;
		 my $mot1=$sa1;
 		  for $j (0..$#plik){
 		   if($kod[$j]=~/^$mot1/){
 		   $sum1=$sum1+$ilosc[$j]; 
 		   }
		  }
			if($sum1 == 0){
			$p1=0;}
			elsif($sum1>0){
			$p1=$sum1/$sum0;
 		   	printf CTOUT "%s\t%7.3f\n",$mot1,$p1;}#print  
		 } 	 		 
		
		close (CTOUT);
		}#nastepny istniejacy sekwencyjny
	  }#nastepny motyw sekwencyjny			
	 }
	}
 }#nastepny katalog z motywami sekwencyjnymi
chdir "$rootdir" or die "Nie moge wejsc do d12: $rootdir $!\n";	
 }#end of subroutine
1;
	
