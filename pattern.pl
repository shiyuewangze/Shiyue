#!/usr/bin/perl -w
#Author: yangxin
#Date:
#Modified:
#Description:Pval/FDR calculated in the levle of whole genome;
use strict;

my $intersectBed="/data/software/bedtools2-2.25.0/bin/intersectBed";
my $peak_file="./Sample_overlap.bed";
my $r="./ratio-position_Homo";
my $intron="./intron-num";
my $utr5="./utr5-num";
my $utr3="./utr3-num";
my $cds="./cds-num";
my $plot="./peak-inter-for-plot";
my $fate=2;
my $LYL=1;
my $label;

system("$intersectBed -a $peak_file -b $r -f 0.51 -wa -wb >$plot");
&m6A_distribution($plot,$r,$fate-1,$cds,$intron,$utr5,$utr3);
for (my $i=0;$i<=$LYL ;$i++) {
	my $j=$i*0.5;
	$j=sprintf("%.1f",$j);
	$label.="\"$j\",";
}
chop $label;
my $limy=sprintf("%.1f",$LYL/2);

my $r_temp ="./plot_tempfileR\.txt";
my $pdf="./Plot_Fig.pdf";
open(R, ">$r_temp");
print R "pdf\(\"$pdf\"\)\n";
#print R "layout(matrix(1:4,2,2))\n";

#print R "par(mar=c(5,4,4,2))\n";
print R "a1<-read.table(\"$utr5\")\n";
print R "b1<-read.table(\"$cds\")\n";
print R "c1<-read.table(\"$utr3\")\n";
print R "x1<-c(a1\$V1,b1\$V1+100,c1\$V1+200)\n";
print R "y1<-c(a1\$V2,b1\$V2,c1\$V2)\n";
print R "plot(smooth.spline(x1,y1,df=50)\$x,smooth.spline(x1,y1,df=50)\$y,type=\"l\",xaxt=\"n\",yaxt=\"n\",xlab=\"\",ylab=\"\",bty=\"n\",ylim=c(0,1.2))\n";
print R "mtext(\"5'UTR\",side=1,line=1.2,at=50,font=2,cex=0.9)\n";
print R "mtext(\"CDS\",side=1,line=1.2,at=150,font=2,cex=0.9)\n";
print R "mtext(\"3'UTR\",side=1,line=1.2,at=250,font=2,cex=0.9)\n";#font=2:fold, font=4:fold&italy
print R "mtext(expression(paste(bold(Percentage),\" \",bold(of),\" \",bold(m^\"6\"),bold(A),\" \",bold(Peaks),\" \",bold(\"(%)\"))),side=2,las=0,line=2.5,cex=0.9)\n";
print R "axis(1,at=c(0,99,199,299),label=c(\"\",\"\",\"\",\"\"))\n";
print R "axis(2,at=c(0:$LYL)*0.5,label=c($label),las=1)\n";
print R "abline(v=99,lty=2,col=\"red\")\n";
print R "abline(v=199,lty=2,col=\"red\")\n";
print R "mtext(\"B\",side=3,font=2,line=1,at=-90)\n";

print R "dev\.off\(\)\n";
close R;
system ("R --vanilla --file=\"$r_temp\"");
sub max{
	my ($a,$b)=@_;
	if ($a>=$b) {
		return $a;
		}
	else{
		return $b;
		}
	}

sub min{
	my ($a,$b)=@_;
	if ($a<=$b) {
		return $a;
	}
	else{
		return $b;
	}
}



sub m6A_distribution {
	my ($pifp,$rp,$fra,$cds,$intron,$utr5,$utr3)=@_;
	open IN,"< $pifp";
	my %hash;
	while (my $aline=<IN>) {
		chomp $aline;
		my @peak=split /\t/,$aline;
		if (@peak>4) {
			push @{$hash{$peak[7]}},$aline;
		}
	}
	close IN;
	open IN2,"< $rp";
	my (%hc,%hi,%hu5,%hu3);
	my $total=0;
	my ($pc,$pi,$pu5,$pu3)=(0,0,0,0);
	#print "*************************\n";
	while (my $a=<IN2>) {
		chomp $a;
		my @pos=split /\t/,$a;
		if (defined $hash{$pos[3]}) {
			foreach my $i (0..$#{$hash{$pos[3]}}) {
				my @arr=split /\t/,$hash{$pos[3]}[$i];
				my $m=&max($arr[1],$pos[1]);
				my $n=&min($arr[2],$pos[2]);
				my ($rs,$re);
				
				if ($pos[-1]=~/\+/) {
					$rs=int((($m-$pos[1])/$pos[4]+$pos[5])*100);
					$re=int((($n-$pos[1])/$pos[4]+$pos[5])*100);
				}
				if ($pos[-1]=~/\-/) {
					$rs=int((($pos[2]-$n)/$pos[4]+$pos[5])*100);
					$re=int((($pos[2]-$m)/$pos[4]+$pos[5])*100);
				}
				#print "$pos[3]\t$pos[1]\t$pos[2]\t$pos[5]\t$pos[6]\t$arr[1]\t$arr[2]\t$rs\t$re\n";
				if ($pos[3]=~/cds/){
						$total++;
						if ($fra==0){
						$pc=int(($rs+$re)/2);
						$hc{$pc}++;
						
						}
					else{
						for (my $s=0;$s<=$fra;$s++) {
							$pc=$rs+int($s*($re-$rs)/$fra);
							$hc{$pc}++;
						}
					}
				}
				if ($pos[3]=~/intron/){
					
					if ($fra==0){
						$pi=int(($rs+$re)/2);
						$hi{$pi}++;
						}
					else{
						for (my $s=0;$s<=$fra;$s++) {
							$pi=$rs+int($s*($re-$rs)/$fra);
							$hi{$pi}++;
						}
					}
				}
					if ($pos[3]=~/utr5/) {
						$total++;
						if ($fra==0){
						$pu5=int(($rs+$re)/2);
						$hu5{$pu5}++;
						}
					else{
						foreach my $s (0..$fra) {
							$pu5=$rs+int($s*($re-$rs)/$fra);
							$hu5{$pu5}++;
							}
						}
					}
					if ($pos[3]=~/utr3/) {
						$total++;
						if ($fra==0){
							$pu3=int(($rs+$re)/2);
								$hu3{$pu3}++;
					}
					else{
						foreach my $s (0..$fra) {
							$pu3=$rs+int($s*($re-$rs)/$fra);
							$hu3{$pu3}++;
						}
					}
				}
			}
		}
	}
	close IN2;
	$total*=($fra+1);
	#print "total:$total\n";
	my (@cds,@intron,@utr5,@utr3);
	
	open CDS,"> $cds";
	delete $hc{0},if (exists $hc{0});
	delete $hi{0},if (exists $hi{0});
	delete $hu5{0},if (exists $hu5{0});
	delete $hu3{0},if (exists $hu3{0});
	foreach my $key(sort{$a<=>$b}keys %hc) {
			my $c=$hc{$key}/$total*100;
		printf CDS "$key\t%.2f\n",$c;
		#print "$key\t$hc{$key}\t$c\n";
		push @cds,$c;
	}
	close CDS;
	open INTRON,"> $intron";
	foreach my $key(sort{$a<=>$b}keys %hi) {
			my $c=$hi{$key}/$total*100;
		printf INTRON "$key\t%.2f\n",$c;
		#print "$key\t$hc{$key}\t$c\n";
		push @intron,$c;
	}
	close INTRON;
	open UTRF,"> $utr5";
	foreach my $key(sort{$a<=>$b}keys %hu5) {
			my $c=$hu5{$key}/$total*100;
		printf UTRF "$key\t%.2f\n",$c;
		#print "$key\t$hu5{$key}\t$c\n";
		push @utr5,$c;
	}
	close UTRF;
	open UTRT,"> $utr3";
	foreach my $key(sort{$a<=>$b}keys %hu3) {
			my $c=$hu3{$key}/$total*100;
		printf UTRT "$key\t%.2f\n",$c;
		#print "$key\t$hu3{$key}\t$c\n";
		push @utr3,$c;
	}
	close UTRT;
	@cds=sort {$a<=>$b} @cds;
	@intron=sort {$a<=>$b} @intron;
	@utr5=sort {$a<=>$b} @utr5;
	@utr3=sort {$a<=>$b} @utr3;
	my $m1=&max($intron[-1],$utr5[-1]);
	my $m2=&max($m1,$utr3[-1]);
	my $m3=&max($m2,$cds[-1]);
	$LYL=int($m3/0.5)+1;
	#print"$intron[-1]\t$utr5[-1]\t$utr3[-1]\t$LYL\n";
}
