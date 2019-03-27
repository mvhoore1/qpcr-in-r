# qpcr-in-r
A script to analyse qPCR results using R and the EasyqpcR package



What we want the script to do, assuming the input file is the one in the mock data. 

# Transform ct_vals by doing 2^-X
# this "2" should be able to change at some point because different primers will have different effciencies. However this is only important when you compare different primers/genes with each other in absolute terms

#group the biological replicates of the reference genes, maybe in a separate table??

# divide every sample with the, just grouped reference genes whose xp_factor and biol_rep match.
#this is to control for the amount of RNA in each sample (reference gene is assumed to be stable)

#now it needs to be normalized to the WT.  Group the wild-type per gene. 
#divide every sample by the just grouped wt of the same gene.
# at this point all the reference gene samples should, when averaged, be 1. And all the WT samples of every gene when averaged should be 1



#Now here the program can branch. Because we want to either
A see what all the different biol_reps do
B see what all the xp_factors do

# if  A. get mean for every biol_rep and get STDEV for every Biol_rep. And then be able to make a graph of al the Biol_reps of a gene with their STDEV
# if  B. get mean for every xp_factor and get STDEV for every xp_factor . And then be able to make a graph of al the xp_factors of a gene with their STDEV
