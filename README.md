# acc2tree
This script is under development and more features are going to be added.

It takes user's accession ID and utilize a blastp search in NCBI database using the user's defined maximum target sequences and makes a list of chosen accession IDs under these conditions:
1- Ids with lowest e-value
2- If the ref seq. is available chooses that along the above condition.
3- If the organism has only one accession ID, gets that no matter what its e-value is.
The script then uses the FASTA file made from the chosen sequences to perform alignment using MUSCLE.
The aligned sequences file will be used to make a phylogenetic tree with IQtree or RAxML based on the user's choice in the begining of the process.
The MUSCLE, IQtree and RAxML should be installed on the system and their directory should be added to the system PATH.
The default options for the tree making programs are as follows:
RAxML: "-m", "PROTGAMMAAUTO", "-f", "a", "-x", "1000", "-p", "1988", "-k", "-N", "autoMRE"
IQtree:"-m", "MFP", "-nt", "AUTO", "-bb", "1000", "-alrt", "1000", "-quiet"
