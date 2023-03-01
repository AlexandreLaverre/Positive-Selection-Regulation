use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readFasta{

    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];
    
    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input,$path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script checks if we have the right genome assembly for a MAF file.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

my %parameters;

$parameters{"species"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathMAF"}="NA";
$parameters{"minAlnSize"}="NA";
$parameters{"maxNbAln"}="NA";
$parameters{"pathOutput"}="NA";
my @defaultpars=("species","pathGenomeSequence","pathMAF", "minAlnSize", "maxNbAln", "pathOutput");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

##############################################################
##############################################################

print "Reading genome sequence...\n";
my %genome;

readFasta($parameters{"pathGenomeSequence"}, \%genome);

print "Done.\n";

##############################################################

print "Extracting chromosome sizes...\n";

my %genomechrsizes;

foreach my $chr (keys %genome){
    my $size=length $genome{$chr};

    if(exists $genomechrsizes{$size}){
	push(@{$genomechrsizes{$size}}, $chr); 
    } else{
	$genomechrsizes{$size}=[$chr];
    }
}

print "Done.\n";

##############################################################

print "Reading MAF and checking sequences...\n";

my $sp=$parameters{"species"};
my $minsize=$parameters{"minAlnSize"}+0;
my $maxnb=$parameters{"maxNbAln"};

my $nbdone=0;

my $input;
my $pathmaf=$parameters{"pathMAF"};
my @s=split("\\.");
my $ext=$s[-1];

if($ext eq "gz"){
    open($input, "zcat $pathmaf |");
} else{
    open($input, $pathmaf);
}

my $line=<$input>;

my %chromonotfound;
my %chromocorresp;
my %nbok;
my %nbnotok;
my %chrsizesaln;

while($line && (($maxnb eq "NA") || ($nbdone<=$maxnb))){
    my $prefix=substr $line, 0, 1;

    if($prefix eq "s"){
	chomp $line;
	$line =~ tr/ //s;
	my @s=split("\t", $line);
	my @spchr=split("\\.", $s[1]);

	if($spchr[0] eq $sp){
	    shift @spchr;
	    my $chr=join(".", @spchr);
	    my $size=$s[3]+0;

	    my $thischrsize=$s[5]+0;
	    	    
	    if($size >= $minsize){
		my $strand=$s[4];
	 
		if($strand eq "+"){
		    $chrsizesaln{$chr}=$thischrsize;
				   
		    my $start=$s[2]+0;

		    my $synchr="NA";

		    my %possiblesyn;

		    if(exists $chromocorresp{$chr}){
			$synchr=$chromocorresp{$chr};
		    } else{
			if(exists $genome{$chr}){
			    $synchr=$chr;
			} else{
			    if(exists $genome{"chr".$chr}){
				$synchr="chr".$chr;
			    } else{
				my $pref=substr $chr, 0, 3;
				if($pref eq "chr"){
				    my $sc=substr $chr, 3;
				    
				    if(exists $genome{$sc}){
					$synchr=$sc;
				    } else{
					## check by size
									
					if(exists $genomechrsizes{$thischrsize}){
					    my $nbsyn=@{$genomechrsizes{$thischrsize}};

					    if($nbsyn==1){
						$synchr=${$genomechrsizes{$thischrsize}}[0];
					    } else{
						foreach my $c (@{$genomechrsizes{$thischrsize}}){
						    $possiblesyn{$c}=1;
						}
					    }
					}
				    }
				} else{
				    if(exists $genomechrsizes{$thischrsize}){
					my $nbsyn=@{$genomechrsizes{$thischrsize}};

					if($nbsyn==1){
					    $synchr=${$genomechrsizes{$thischrsize}}[0];
					} else{
					    foreach my $c (@{$genomechrsizes{$thischrsize}}){
						$possiblesyn{$c}=1;
					    }
					}
				    }
				}
			    }
			}
		    }
		    
		    if($synchr ne "NA"){
			$chromocorresp{$chr}=$synchr;

			if(!exists $nbok{$chr}){
			    $nbok{$chr}=0;
			}

			if(!exists $nbnotok{$chr}){
			    $nbnotok{$chr}=0;
			}
			
			my $sequence=uc $s[6];
			$sequence =~ s/-//gi;

			my $correspseq=uc (substr $genome{$synchr}, $start, $size);

			if($sequence eq $correspseq){
			    $nbok{$chr}++;
			    
			} else{
			    print "Different sequences found for ".$chr.":\n";
			    print $correspseq."\n";
			    print $sequence."\n";
			    
			    $nbnotok{$chr}++;
			}
		    } else{
			my $nbpos=keys %possiblesyn;

			if($nbpos>=1){
			    my $sequence=uc $s[6];
			    $sequence =~ s/-//gi;
			    			    
			    foreach my $sc (keys %possiblesyn){
				
				my $correspseq=uc (substr $genome{$sc}, $start, $size);
				
				if($sequence eq $correspseq){
				    $chromocorresp{$chr}=$sc;
				    
				    if(!exists $nbok{$chr}){
					$nbok{$chr}=1;
				    }
				    
				    if(!exists $nbnotok{$chr}){
					$nbnotok{$chr}=0;
				    }

				    last;
				}
			    }
			} else{
			    print "Cannot find correspondence for chr ".$chr." size ".$thischrsize."\n";
			}
		    }
		    
		    $nbdone++;

		    if($nbdone%10000==0){
			print $nbdone." alignments done.\n";
		    }
		}
	    }
	}
    }

     $line=<$input>;
}

close($input);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "ChrMAF\tChrMAFSize\tChrAssembly\tChrAssemblySize\tNbOK\tNbKO\n";

foreach my $mafchr (keys %chrsizesaln){
    my $mafsize=$chrsizesaln{$mafchr};

    if(exists $chromocorresp{$mafchr}){
	my $asschr=$chromocorresp{$mafchr};

	my $asssize=length $genome{$asschr};

	print $output $mafchr."\t".$mafsize."\t".$asschr."\t".$asssize."\t".$nbok{$mafchr}."\t".$nbnotok{$mafchr}."\n";
    } else{
	print $output $mafchr."\t".$mafsize."\tNA\tNA\tNA\tNA\n";
    }
}
close($output);

print "Done.\n";

##############################################################
