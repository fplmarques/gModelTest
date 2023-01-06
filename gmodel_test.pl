#!/usr/bin/perl
#
# This program uses Garli to select substitution models for phylogenetic analyeses usin Maximum Likelihood 
# as optimality criterium. As such it requires [Garli](https://code.google.com/archive/p/garli/)
#
#  Usage:
#       $ gmodel_test.pl [option] aligned_sequences.fas
# 
system clear;
use Getopt::Long;
use POSIX qw(ceil);
#use Parallel::ForkManager;
#$pm = new Parallel::ForkManager(8);
############### CONTROLERS
#
#
$start_run = time();
$overall_best_AICc_score = 100000000;
#
#
$fixed_ML_GTR_topology = 0;									# used 1 for fixed topology and 0 for ML optimized for each model
GetOptions ('fixed|fix|f' => \$fixed_ML_GTR_topology);   					# if option "-fix"is imposed test consider fixed topology
$garli_replicates = 3;										# controls number of replicates in Garli during start tree search
$garli_output_prefix_fixed_tree = "starting_fixed_tree";
$garli_output_prefix_ml_optimized = "ml_optimized";
$input_file = "$ARGV[0]";
$input_file =~ s/\.fas//g;
$results_output_file = "$input_file\_model_test_results";
#
#
if ($fixed_ML_GTR_topology =~ /1/){
	$results_output_file = "$results_output_file\_fixed.txt";
	} else {$results_output_file = "$results_output_file\_optimized.txt";}
#
#
open (FILE, "$ARGV[0]");
$Ntaxa = 0;
while (<FILE>) {
	    if ($_ =~ m/^>/){
		$Ntaxa++;
		}
		if ($_ =~ m/^[-ACGTN]/){
		chomp($_);
		$sequence_init = $_;
		$sequence_length = $sequence_length+length($sequence_init);
		}
}
$sequence_length = $sequence_length/$Ntaxa;		# This value is the first estimated sequence length, might contain error.
close (FILE);
#
print "\n\n\n";
print "     ######################################################################\n";
print "     #                         GModelTest v.00                            #\n";
print "     ######################################################################\n\n\n";
print "     ----> Parsing files $ARGV[0] ...\n";
print "     ----> File $ARGV[0] has $Ntaxa terminals.\n";
print "     ----> Sequences are $sequence_length base pairs long.\n\n";
if ($fixed_ML_GTR_topology =~ /1/){
print "     ----> Base tree for likelihood calculations: fixed (option -f imposed) \n";
print "     ----> Fixed topology estimated by GTR+I+G in Garli \n\n";
} else {
print "     ----> Base tree for likelihood calculations = ML tree optimized \n";
print "     ----> To impose fixed topology use option -f \n";
print "     ----> ML tree optimized tests take longer! \n\n";
}
#
########### LOG FILE HEADINGS ###########
open(LOG,">$results_output_file");
print LOG "MODEL \t -lnL\t K \t AIC \t nK ratio \t AICc \t AICc1 \n";
close (LOG);
#########################################
if ($fixed_ML_GTR_topology =~ /1/){
$tree = 0;								# no topology will be estimated during model selection. K does not account for that
print "     ######################################################################\n";
print "     #                       Getting Tree in GARLI                        #\n";
print "     ######################################################################\n\n\n";
print "   ------> Building configuration file for initial run ...\n";
############# writes configuration file for initial run
#
open(OUT,">garli_starting_fixed_tree.conf");
print OUT "# this is a model test file
[general]
datafname = $ARGV[0]
ofprefix = $garli_output_prefix_fixed_tree
outputphyliptree = 1
streefname = stepwise
attachmentspertaxon = 10
constraintfile = none
searchreps = $garli_replicates
outgroup = 1

outputeachbettertopology = 0
outputcurrentbesttopology = 0

enforcetermconditions = 1
genthreshfortopoterm = 5000
scorethreshforterm = 0.001
significanttopochange = 0.01

writecheckpoints = 0
restart = 0

randseed = -1
availablememory = 512
logevery = 10
saveevery = 100
refinestart = 1
outputphyliptree = 0
outputmostlyuselessfiles = 0
collapsebranches = 1

linkmodels = 0
subsetspecificrates = 1

[model0]
#GTR+I+G
datatype = nucleotide
ratematrix = 6rate
statefrequencies = estimate
ratehetmodel = gamma
numratecats = 4
invariantsites = estimate

[master]
bootstrapreps = 0

nindivs = 4
holdover = 1
selectionintensity = 0.5
holdoverpenalty = 0
stopgen = 5000000
stoptime = 5000000

startoptprec = 0.5
minoptprec = 0.01
numberofprecreductions = 2
treerejectionthreshold = 20.0
topoweight = 0.01
modweight = 0.002
brlenweight = 0.002
randnniweight = 0.1
randsprweight = 0.3
limsprweight =  0.6
intervallength = 100
intervalstostore = 5

limsprrange = 6
meanbrlenmuts = 5
gammashapebrlen = 1000
gammashapemodel = 1000
uniqueswapbias = 0.1
distanceswapbias = 1.0


resampleproportion = 1.0
inferinternalstateprobs = 0
\n";
close (OUT);
############# end of wrting garli configuration files
print "   ------> Running GARLI to get tree. Garli will be doing $garli_replicates replicate searches, it could take some time  ...\n";
############## run garli
system "garli garli_starting_fixed_tree.conf  >\/dev\/null 2>\/dev\/null";
#print "I am trying to run garli"; 
############## end of garli runs
print "   ------> Got the tree! Now let me use it!\n\n\n";
#
############# Getting number of terminais for K calculation
#open (START_LOG, "$garli_output_prefix_fixed_tree.screen.log");
#while (<START_LOG>) {
#	if ($_ =~ m/^\s\s+(\d+)\ssequences./){
#		$Ntaxa=$1;
#	}
#}
#close (START_LOG);
print "the dataset has $Ntaxa terminals.";
######################################## HERE WE FINISH THE FIRST RUN 
} else {$tree = 1;} # this is the conditional operator for fixed topologies, also controlls values of K
#
print "     ######################################################################\n";
print "     #                      STARTING MODEL EVALUATIONS                    #\n";
print "     ######################################################################\n\n\n\n";
print "   ------> There we go  ...\n";

for $ratematrix ("\(a a a a a a\)","\(a b a a b a\)","\(a b a a c a\)","\(a b c c b a\)","\(a b a c b c\)","\(a b c a b c\)","\(a b c c d a\)","\(a b a c d c\)","\(a b c a d c\)","\(a b c d b e\)","\(a b c d e f\)"){
	#my $pid = $pm->start and next;
	for $statefrequencies ("equal","estimate"){
		for $invariantsites ("none","estimate"){
			for $ratehetmodel ("none","gamma"){
				if ($ratehetmodel =~ /gamma/){
						$numratecats = 4;
				} else {$numratecats = 1;}				
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "JC";
					$free_parameters=0;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "JC+I";
					$free_parameters=1;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "JC+G";
					$free_parameters=1;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "JC+I+G";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "F81";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "F81+I";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "F81+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a a a a a a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "F81+I+G";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "K80";
					$free_parameters=1;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "K80+I";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "K80+G";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "K80+I+G";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "HKY";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "HKY+I";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "HKY+G";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b a a b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "HKY+I+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TrNef";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TrNef+I";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TrNef+G";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TrNef+I+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TrN";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TrN+I";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TrN+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b a a c a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TrN+I+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM1";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM1+I";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM1+G";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM1+I+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM1uf";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM1uf+I";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM1uf+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c c b a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM1uf+I+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM2";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM2+I";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM2+G";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM2+I+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM2uf";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM2uf+I";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM2uf+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b a c b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM2uf+I+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM3";
					$free_parameters=2;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM3+I";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM3+G";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM3+I+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM3uf";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TPM3uf+I";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM3uf+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c a b c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TPM3uf+I+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM1ef";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM1ef+I";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM1ef+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM1ef+I+G";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM1";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM1+I";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM1+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c c d a\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM1+I+G";
					$free_parameters=8;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM2ef";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM2ef+I";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM2ef+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM2ef+I+G";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM2";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM2+I";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM2+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b a c d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM2+I+G";
					$free_parameters=8;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM3ef";
					$free_parameters=3;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM3ef+I";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM3ef+G";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM3ef+I+G";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM3";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TIM3+I";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM3+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c a d c\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TIM3+I+G";
					$free_parameters=8;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TVMef";
					$free_parameters=4;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TVMef+I";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TVMef+G";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TVMef+I+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "TVM";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "TVM+I";
					$free_parameters=8;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TVM+G";
					$free_parameters=8;
				}
				if (($ratematrix =~ /\(a b c d b e\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "TVM+I+G";
					$free_parameters=9;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "SYM";
					$free_parameters=5;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "SYM+I";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "SYM+G";
					$free_parameters=6;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /equal/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "SYM+I+G";
					$free_parameters=7;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /none/)){
					$model_name = "GTR";
					$free_parameters=8;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /none/)){
					$model_name = "GTR+I";
					$free_parameters=9;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /none/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "GTR+G";
					$free_parameters=9;
				}
				if (($ratematrix =~ /\(a b c d e f\)/) && ($statefrequencies =~ /estimate/) && ($invariantsites =~ /estimate/) && ($ratehetmodel =~ /gamma/)){
					$model_name = "GTR+I+G";
					$free_parameters=10;
				}
###### CALCULAING K
#
# Optimized free parameters (K) = substitution parameters + N branch lengths + topology
# Branches are a function of 2n-3 in which n is the number of terminals
#
$NumberOfParameters_K = ($free_parameters+((2*$Ntaxa)-3)+$tree);
$scores{$model_name}{'NumberOfParameters_K'} = $NumberOfParameters_K;
#
#
############# after selecting a model writes configuration file
#
## Heading of configuration file depends on choice of trees for model selection
# for fixed topologies
#
if ($fixed_ML_GTR_topology =~ /1/){
open(OUT,">fixed_tree_model_tests.conf");
print OUT "# this is a model test file
[general]
datafname =  $ARGV[0]
streefname = $garli_output_prefix_fixed_tree.best.tre
ofprefix = $garli_output_prefix_fixed_tree
optimizeinputonly = 1
streefname = stepwise
attachmentspertaxon = 10
outgroup = 1
\n";
close (OUT);
} # conditional for heading
#
# for Ml optimized
#
if ($fixed_ML_GTR_topology =~ /0/){
open(OUT,">ml_optimized_model_test.conf");
print OUT "# this is a model test file
[general]
datafname = $ARGV[0]
ofprefix = $garli_output_prefix_ml_optimized
streefname = stepwise
attachmentspertaxon = 10
constraintfile = none
searchreps = $garli_replicates
outgroup = 1
\n";
close (OUT);
}
#
# controls second part of Garli configuration file
#
if ($fixed_ML_GTR_topology =~ /1/){
open(OUT,">>fixed_tree_model_tests.conf");
} else {open(OUT,">>ml_optimized_model_test.conf");}
print OUT "
outputeachbettertopology = 0
outputcurrentbesttopology = 0

enforcetermconditions = 1
genthreshfortopoterm = 5000
scorethreshforterm = 0.001
significanttopochange = 0.01

writecheckpoints = 0
restart = 0

randseed = -1
availablememory = 512
logevery = 10
saveevery = 100
refinestart = 1
outputphyliptree = 0
outputmostlyuselessfiles = 0
collapsebranches = 1

linkmodels = 0
subsetspecificrates = 1

[model0]
#$model_name
datatype = nucleotide
ratematrix = $ratematrix
statefrequencies = $statefrequencies
ratehetmodel = $ratehetmodel
numratecats = $numratecats
invariantsites = $invariantsites

[master]
bootstrapreps = 0

nindivs = 4
holdover = 1
selectionintensity = 0.5
holdoverpenalty = 0
stopgen = 5000000
stoptime = 5000000

startoptprec = 0.5
minoptprec = 0.01
numberofprecreductions = 2
treerejectionthreshold = 20.0
topoweight = 0.01
modweight = 0.002
brlenweight = 0.002
randnniweight = 0.1
randsprweight = 0.3
limsprweight =  0.6
intervallength = 100
intervalstostore = 5

limsprrange = 6
meanbrlenmuts = 5
gammashapebrlen = 1000
gammashapemodel = 1000
uniqueswapbias = 0.1
distanceswapbias = 1.0


resampleproportion = 1.0
inferinternalstateprobs = 0
\n";
close (OUT);
############# end of writting garli configuration files
print "   ------> Configuration file for $model_name has been written, I will show you:\n";
print "\t\t\t Model Name = $model_name\n";
print "\t\t\t Rate Matrix = $ratematrix\n";
print "\t\t\t State Frequencies = $statefrequencies\n";
print "\t\t\t Invariant Sites = $invariantsites\n";
print "\t\t\t Rate Het. Model = $ratehetmodel\n";
print "\t\t\t Num. Rate Cats. = $numratecats\n";
print "\t\t\t Num. Free. Parameters = $free_parameters\n";
print "\t\t\t K = $NumberOfParameters_K\n\n";
############## run garli on models
if ($fixed_ML_GTR_topology =~ /1/){
$garli_start_run = time();
print "     ----> Optimizing fixed tree using GARLI under $model_name model ...\n";
system "garli fixed_tree_model_tests.conf >\/dev\/null 2>\/dev\/null";
} # conditional for getting scores under fixed topology
if ($fixed_ML_GTR_topology =~ /0/){
$garli_start_run = time();
print "     ----> Searching best fit in GARLI under $model_name model ...\n";
system "garli ml_optimized_model_test.conf >\/dev\/null 2>\/dev\/null";
} # conditional for getting scores under fixed topology
#
#
######################### end of garli runs
######################### getting best score from Garli for a given model ##################
## for fixed trees
#
if ($fixed_ML_GTR_topology =~ /1/){
open (FILE, "$garli_output_prefix_fixed_tree.screen.log");
while (<FILE>) {
	if ($_ =~ m/^Final score \= (\-\d+\.\d+)\n/){
		$lnL_score=$1;
	}
}
close (FILE);
$Best_nlL_Score_found = $lnL_score;
$scores{$model_name} = $Best_nlL_Score_found;
} # conditional for getting scores under fixed topology
#
## for ML optimized tress
if ($fixed_ML_GTR_topology =~ /0/){
open (FILE, "$garli_output_prefix_ml_optimized.log00.log");
while (<FILE>) {
	if ($_ =~ m/^Final\s*(\-\d+\.\d+).*/){
		$lnL_score=$1;
		$count_replicates++;
		print "                   Score for replicate $count_replicates: $lnL_score \n";
		push(@lnL_scores_found, $lnL_score=$1)
	}
}
close (FILE);
####################################################################
open (FILE, "ml_optimized.screen.log");
while (<FILE>) {
	if ($_ =~ m/^\s+(\d+)\stotal characters.\n/){
		$sequence_length=$1;
		print "                   Total number of characters (sequence length): $sequence_length \n";
	}
	if ($_ =~ m/^\s+(\d+)\sunique patterns in compressed data matrix.\n/){
		$unique_patterns=$1;
		print "                   Unique patterns in compressed data matrix: $unique_patterns \n";
	}
}
close (FILE);
######################################################

@lnL_scores_found = sort(@lnL_scores_found);
$Best_nlL_Score_found = @lnL_scores_found[0];
$scores{$model_name} = $Best_nlL_Score_found;
@lnL_scores_found = ();
$count_replicates = 0;
} # conditional for ML optimized calculations
#
## calculating AIC (AIC = −2l + 2K) and printing -nlL 
#
$AIC_score = ((-2*$Best_nlL_Score_found)+(2*$NumberOfParameters_K));
$AICc_denominator = ($sequence_length-$NumberOfParameters_K-1);		# this is to avoid illegal division by zero
if ($AICc_denominator =~ /0/){
   $AICc_denominator = 1;
}
$AICc_score = $AIC_score+(((2*$NumberOfParameters_K)*($NumberOfParameters_K-1))/$AICc_denominator);
$AICc1_denominator = ($unique_patterns-$NumberOfParameters_K-1);         # this is to avoid illegal division by zero
if ($AICc1_denominator =~ /0/){
   $AICc1_denominator = 1;
}
$AICc1_score = $AIC_score+(((2*$NumberOfParameters_K)*($NumberOfParameters_K-1))/$AICc1_denominator);
$nKratio = $sequence_length/$NumberOfParameters_K;
$scores{$model_name}{'NumberOfParameters_K'} = $NumberOfParameters_K;
$scores{$model_name}{'AIC_score'} = $AIC_score;
$scores{$model_name}{'AICc_score'} = $AICc_score;
$scores{$model_name}{'AICc1_score'} = $AICc1_score;
$scores{$model_name}{'nKratio'} = $nKratio;
#
# reporting
$progress++;
############### calculating total run time
# for garli
$garli_end_run = time();
$garli_run_time = $garli_end_run - $garli_start_run;
$garli_run_time_min = ($garli_run_time/60)%60;
$garli_run_time_sec = $garli_run_time%60;
# running time so far
$partial_run = time();
$run_time_so_far = $partial_run - $start_run;
$run_time_so_far_d = int($run_time_so_far/(24*60*60));
$run_time_so_far_hrs = ($run_time_so_far/(60*60))%24;
$run_time_so_far_min = ($run_time_so_far/60)%60;
$run_time_so_far_sec = $run_time_so_far%60;
# estimated time to complete
$sum_of_garli_run_time = $sum_of_garli_run_time + $garli_run_time;
$average_garli_run_time = $sum_of_garli_run_time/$progress;
$number_of_models_left = 88 - $progress;
$expected_time_to_complition = $average_garli_run_time * $number_of_models_left;
$expected_time_to_complition_d = int($expected_time_to_complition/(24*60*60));
$expected_time_to_complition_hrs = ($expected_time_to_complition/(60*60))%24;
$expected_time_to_complition_min = ($expected_time_to_complition/60)%60;
$expected_time_to_complition_sec = $expected_time_to_complition%60;
############## end of calculating total run time
print "                   Total number of free parameters (including gbranch lengths): $NumberOfParameters_K \n";
print "       --> Best score for this model: $Best_nlL_Score_found \n";
print "       --> Akaike Information Criterion: $AIC_score \n";
print "       --> Corrected Akaike Information Criterion: $AICc_score \n";
print "       --> Corrected Akaike Information Criterion (alt): $AICc1_score \n";
print "       --> Ratio of sample size (n) compared to the number of parameters: $nKratio \n";
print "       --> Computational time for the model: $garli_run_time_min min $garli_run_time_sec s\n";
print "       --> Accumulated running time: $run_time_so_far_d d: $run_time_so_far_hrs hrs: $run_time_so_far_min min: $run_time_so_far_sec secs \n";
print "       --> Expected time to complete accessing all models: $expected_time_to_complition_d d: $expected_time_to_complition_hrs hrs: $expected_time_to_complition_min min: $expected_time_to_complition_sec secs \n\n";
#############################  printing progress report
$gone = ceil(1.136363636*$progress);
print "       --> Progress: |";
for $i (1..$progress){
   print "|";
}
for  $i (1..88-$progress){
   print " ";
}
print "| $gone\/100%\n\n";
#############################
open(LOG,">>$results_output_file");
print LOG "$model_name\t$Best_nlL_Score_found\t$NumberOfParameters_K\t$AIC_score\t$nKratio\t$AICc_score\t$AICc1_score\n";
close (LOG);
##### getting best model
if ("$AICc_score" < "$overall_best_AICc_score"){
	$overall_best_AICc_score = $AICc_score;
	$Name_of_best_model = $model_name;
}
print "       --> The substitution model $Name_of_best_model holds the best AICc so far: $overall_best_AICc_score\n\n";
			}
		}
	}
# $pm->finish;
}
print "     ######################################################################\n";
print "     #                                DONE!                               #\n";
print "     ######################################################################\n\n";
$end_run = time();
$run_time = $end_run - $start_run;
$run_time_d = int($run_time/(24*60*60));
$run_time_hrs = ($run_time/(60*60))%24;
$run_time_min = ($run_time/60)%60;
$run_time_sec = $run_time%60;
print "     ######################################################################\n";
print "     #                         SUMMARY RESULTS                            #\n";
print "     ######################################################################\n\n";
print "       --> Best fitting model (AICc): $Name_of_best_model \n";
print "       --> nlL Score for the model: $scores{$Name_of_best_model} \n";
print "       --> AIC for the model: $scores{$Name_of_best_model}{'AIC_score'} \n"; 
print "       --> AICc for the model: $scores{$Name_of_best_model}{'AICc_score'} \n";
print "       --> gModelTest evaluated 88 models in $run_time_d d: $run_time_hrs h: $run_time_min min: $run_time_sec secs \n\n";
print "       *** Results for all models tested are in $results_output_file \n";
print "       *** Alternative Corrected Akaike Information Criterion (alt) calculations use number of unique patterns in compressed data matrix as sample size.\n";
system `rm -f starting_fixed_tree.*.* ml_optimized*.* garli_starting_fixed_tree.conf fixed_tree_model_tests.conf`;
##### printing log file
open(LOG2,">>GModelTest.log");
print LOG2 "Input File \t Model \t nlL \t K \t AIC \t AICc \t nK Ratio \n";
print LOG2 "$input_file\.fas\t$Name_of_best_model\t$scores{$Name_of_best_model}\t$scores{$Name_of_best_model}{'NumberOfParameters_K'}\t$scores{$Name_of_best_model}{'AIC_score'}\t$scores{$Name_of_best_model}{'AICc_score'}\t$scores{$Name_of_best_model}{'nKratio'}\n";
close (LOG2);
exit();
