
use strict;
#

my %hash=();
my $normalFlag=0;

open(RF,"risk.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($arr[$#arr] eq "low"){
	  $hash{$arr[0]}=0;
	}
	elsif($arr[$#arr] eq "high"){
	  $hash{$arr[0]}=1;
	}
}
close(RF);

my @indexs=();
my @riskArr=();
open(RF,"geneMatrix.txt") or die $!;
open(WF,">sampleSymbol.txt") or die $!;
my $normalCount=0;
my $tumorCount=0;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		print WF "ID";
		for(my $i=1;$i<=$#arr;$i++){
			my @samples=split(/\-/,$arr[$i]);
			my $sampleName=$arr[$i];
			if(exists $hash{$sampleName}){
					push(@indexs,$i);
					push(@riskArr,$hash{$sampleName});
					print WF "\t$arr[$i]";
					$tumorCount++;
					delete($hash{$sampleName});
			}
		}
		print WF "\n";
	}
	else{
		print WF $arr[0];
		foreach my $col(@indexs){
			print WF "\t$arr[$col]";
		}
		print WF "\n";
	}
}
print WF "Risk\t" . join("\t",@riskArr) . "\n";
close(WF);
close(RF);


my $colNum=0;
my $rowNum=0;
%hash=();
my $geneName="Risk";
@indexs=();
my @geneArr=();

open(RF,"sampleSymbol.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		for(my $i=1;$i<=$#arr;$i++){
			my @samples=split(/\-/,$arr[$i]);
			push(@indexs,$i);
			my $sampleName=$arr[$i];
			$hash{$sampleName}=1;
			$colNum++;
		}
	}
	else{
			$rowNum++;
		  if($arr[0] eq $geneName){
			  foreach my $col(@indexs){
				  push(@geneArr,$arr[$col]);
			  }
		  }
	}
}
close(RF);

my $firstGeneVal=$geneArr[0];
my $geneMed=average(@geneArr);

open(RF,"sampleSymbol.txt") or die $!;
open(WF,">$geneName.gct") or die $!;
print WF "#1.2\n";
print WF "$rowNum\t$colNum\n";
open(CLS,">$geneName.cls") or die $!;
print CLS "$colNum\t2\t1\n";
if($firstGeneVal>$geneMed){
	print CLS "#\th\tl\n";
}
else{
	print CLS "#\tl\th\n";
}

#my @samp1e=(localtime(time));
@indexs=();
my @typeArr=();
while(my $line=<RF>){
	chomp($line);

	my @arr=split(/\t/,$line);
	if($.==1){
		print WF "NAME\tDESCRIPTION";
		#if($samp1e[4]>13){next;}
		for(my $i=1;$i<=$#arr;$i++){
			my @samples=split(/\-/,$arr[$i]);
			my $sampleName=$arr[$i];
			if(exists $hash{$sampleName}){
				  push(@indexs,$i);
				  print WF "\t$arr[$i]";
				  #delete($hash{$sampleName});
			  }
		}
		print WF "\n";
	}
	else{
		my $symbolName=$arr[0];
		$symbolName=~s/(.+?)\|.+/$1/g;
		print WF "$symbolName\tna";
		foreach my $col(@indexs){
			print WF "\t$arr[$col]";
		}
			print WF "\n";
		if($arr[0] eq $geneName){
			foreach my $col(@indexs){
				if($arr[$col]>$geneMed){
				  push(@typeArr,"h");
			  }
			  else{
			  	push(@typeArr,"l");
			  }
			}
		}
	}
}
print CLS join("\t",@typeArr) . "\n";
close(WF);
close(CLS);
close(RF);
unlink<"sampleSymbol.txt">;

sub median{  
    my (@data) = sort {$a <=> $b} @_;  #??ԭ??????  
    if(scalar (@data) % 2){  
        return ($data [@data / 2]);  #@data/2?ķ??????ȡ?  
    }else {  
        my ($upper ,$lower);  
        $upper = $data[@data / 2];  
        $lower = $data[@data / 2 -1];  
        return (($lower+$upper) / 2);  
    }  
}  

sub average{
	 my $aver = 0;
	 map{$aver += $_} @_;
	 $aver /=($#_+1);
	 return $aver;
}

