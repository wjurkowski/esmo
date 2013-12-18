# CTcalc - anayzes existing CT
# $Id: CTcalc.pm
package CTcalc;

use strict;
use warnings;
our @ISA = qw(Exporter);
our @EXPORT = qw(calSE);

sub calSE {
my($params)=@_;
my $centrop=$$params{'se_cout'};
my $rentrop=$$params{'se_rout'};
open(SECOUT,">$centrop") or die "Can’t write entropies for AA motifs (columns)\n";
open(SEROUT,">$rentrop") or die "Can’t write entropies for SA motifs (rows)\n";	
print SEROUT "row #\tmotif\tmax p\tSE2\tSEmax2\tdSE2\tSE20\tSEmax20\tdSE20\n";
print SECOUT "col #\tmotif\tmax p\tSE2\tSEmax2\tdSE2\tSE7\tSEmax7\tdSE7\n";
my @aatab=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
my @sctab=('A','B','C','D','E','F','G');
my $aalibry="$$params{'libaadir'}";
my $salibry="$$params{'libsadir'}";
my $cmn=0;
my $rootdir=`pwd`;#katalog roboczy
chomp $rootdir;
my ($kat,$AA2,$AA3,$AA4,$i,@samot,@rcnt,$maxpr,$maxprsa,@motsum,@SEC2,@SEC7,@SEC20,@pr,@SECMAX2,@SECMAX7,@dSEC2,@dSEC7,@MPSA);

chdir "$aalibry" or die "Nie moge wejsc do katalogu d1 $aalibry $!\n";
 foreach $kat(@aatab){#motywy sekwencyjne
 	chdir "$kat" or die "Nie moge wejsc do katalogu d2 $kat $!\n";
	foreach $AA2(@aatab){
  	 foreach $AA3(@aatab){
	  foreach $AA4(@aatab){
		my $motyw=$kat.$AA2.$AA3.$AA4;
		if( -f "$motyw"){#jezeli istnieje otworz go
		$cmn++;#numer kolumny
		open(CTAAINP,"< $motyw") or die "Can�t open SA motifs file of AA motif: $!, $motyw\n";
		my @plik=<CTAAINP>;#wczytaj do tablicy
		chomp @plik;
		close (CTAAINP);
		 my $rmn=0;
		 $motsum[$cmn]=0;
		 $SEC2[$cmn]=0;
		 $SEC7[$cmn]=0;
		 for $i (0..$#plik){#strukturalny dla danego sekwencyjnego
		 $rmn++;#liczba stanow dla danej kolumny
		 my @para=split(/\s+/,$plik[$i]);
		 $samot[$i]=$para[0];
		 $rcnt[$i]=$para[1];
		 $motsum[$cmn]=$motsum[$cmn]+$para[1];}#nastepny strukturalny
		 $maxpr=0;
		 for $i (0..$#plik){#prawdopodobienstwa
		 $pr[$i]=$rcnt[$i]/$motsum[$cmn];
			if($maxpr < $pr[$i]){
			$maxpr=$pr[$i];
			$maxprsa=$samot[$i];}
		 my $entr2=$pr[$i]*(log($pr[$i])/log(2));
		 my $entr7=$pr[$i]*(log($pr[$i])/log(7));
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
	
 my $rmn=0;
 my (@aamot,@ccnt,$SA2,$SA3,$SA4,$maxpraa,@SER2,@SER20,@SERMAX2,@SERMAX20,@MPAA,@dSER2,@dSER20);
chdir "$rootdir/$salibry" or die "Nie moge wejsc do katalogu d3 $rootdir/$salibry $!\n";
 foreach $kat(@sctab){#motywy strukturalne
 	chdir "$kat" or die "Nie moge wejsc do katalogu d4 $kat $!\n";
	foreach $SA2(@sctab){
  	 foreach $SA3(@sctab){
	  foreach $SA4(@sctab){
		my $motyw=$kat.$SA2.$SA3.$SA4;
		if( -f "$motyw"){#jezeli istnieje otworz go
		$rmn++;#numer wiersza
		open(CTSAINP,"< $motyw") or die "Can�t open SA motifs file of AA motif: $!, $motyw\n";
		my @plik2=<CTSAINP>;#wczytaj do tablicy
		chomp @plik2;
		close (CTSAINP);
		 my $cmn=0;
		 $motsum[$rmn]=0;
		 $SER2[$rmn]=0;
		 $SER20[$rmn]=0;
		 for $i (0..$#plik2){#sekwencyjny dla danego strukturalnego
		 $cmn++;#liczba stanow dla danego wiersza
		 my @para=split(/\s+/,$plik2[$i]);
		 $aamot[$i]=$para[0];
		 $ccnt[$i]=$para[1];
		 $motsum[$rmn]=$motsum[$rmn]+$para[1];}#nastepny sekwencyjny
		 $maxpr=0;
		 for $i (0..$#plik2){#prawdopodobienstwa
		 $pr[$i]=$ccnt[$i]/$motsum[$rmn];
			if($maxpr < $pr[$i]){
			$maxpr=$pr[$i];
			$maxpraa=$aamot[$i];}
		 my $entr2=$pr[$i]*(log($pr[$i])/log(2));
		 my $entr20=$pr[$i]*(log($pr[$i])/log(20));
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
 	}#nastepny katalog z motywami strukturalnymi
}#end of subroutine