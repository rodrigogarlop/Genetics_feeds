# IMPORTANT NOTE: The script was originally intended for use with absolute values in the input table. It has been changed to use relative values but frequency will be kept as a fraction calculated from the total sum of abundances. Please consider this when selecting thresholds. By doing this, the min_reads is actually the adequate filter to remove by relative abundance
# Update 2021-03-18: Removed a min observations of 2 since it may now be used for relative abundance. Also, groups are now optional and I defaulted to 1 group = nrow(). The proportion of min_samples can now be calculated with 
# Update 2020-06-01: Sample singletons are no longer considered (1s in any part of the table) in order to emulate ASV creation behaviour in dada2. I also added a min_reads filter (check position). This can be used to filter out singletons
# Update 2020-02-25: The script was adapted to use a normal contingency table with NO taxonomy in the last column (it may be the first, though, just for the name of each line). I also removed the check for # contructed from BIOM
# Sparsity_reduction_merged_tables.R
# 2019-09-06 By Rodrigo García-López
# This script was tested with R 3.6.1 (2019-07-05) -- "Action of the Toes"
# It is intended to take a contingency table created from a biom file and filter features (given a % cutoff for rows) according to their minimum frequency and minimum samples (min columns) in order to reduce data sparsity (limit the core features explored).
# It was created because QIIME 1 and 2's scripts are not intended for filtering features from different tables. This arises from the fact that prefiltering each tables permanently removes those features that do not pass the filters so that they are not recoverable afterwards, even if the information from the other tables that are merged effectively "saved" those items. This should be avoided as it creates arbitrary null observations for those features appearing in other tables and may inflate differentially distributed items by eliminating observations where there were actually any.
# The input must be a single otu table (output from converting biom tables) from stdin where columns appear in the order of the tables.The first column has the OTU ids of taxa.
# Additionally, parameters are received as arguments for the min frequency cutoff, the min sample cutoff and the number of items per sample (each separately).
# Run as follows:
# cat table.tsv|Rscript Sparsity_reduction_merged_tables.R <min_freq> <min_samples> <output_file> <#_samples_table_1> [#_samples_table_1] ... [#_samples_table_n]

# Tested with command:
# cat 02_split_tables/genetics_lvl4-H_grp_A.tsv|Rscript Sparsity_reduction_merged_tables.R 1 0.00005 0.25 out.tsv 6 6 6 6
# df <- read.table("02_split_tables/genetics_lvl4-H_grp_A.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names = FALSE)
# min_reads = 0
# min_freq = 0
# min_samples = 0.5
# output_name = "test"
# group_sizes = 6

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<4) { # at least, five arguments are required
  stop("A minimum of 5 arguments is mandatory: cat table.tsv|Rscript Sparsity_reduction_merged_tables.R <min_reads> <min_freq> <min_samples> <output_file> [#_samples_table_1] ... [#_samples_table_n]", call.=FALSE)
}
min_reads <- as.numeric(args[1]) # Get the minimum reads to consider valid (per row)
min_freq <- as.numeric(args[2]) # Get the cutoff for minimum frequence in terms of proportion (0-1; recommended: 0.0001) (this is applyied to the whole table after filtering
min_samples <- as.numeric(args[3]) # Get the minimum proportion of samples presenting the feature (0-1; recommended: 25%)
output_name <- as.character(args[4]) # Get a prefix for the output (can include paths if existent)
group_sizes <- as.numeric(args[5:length(args)]) # Create a vector with the total group sizes (variable size and optional)
groups <- length(group_sizes) # Store the number of total groups to consider
paste("Min freq cutoff:", min_freq, "  Min samples cutoff:", min_samples, "  Group size:",unlist(group_sizes))
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names = FALSE) # set skip to 1 if the "Constructed from biom file" line is present
df <- rowsum(df[2:length(df)], group=df[[1]]) #sum duplicate items by feature (row) name (if any) UNCOMMENT if required
df <- df[rowSums(df)>=min_reads,] # Update 2020-06-01: Remove all 1s from the table and delete resulting empty rows (we can't trust unique sample singletons as they may well arise from errors). This is intended to filter those with at least 2 reads but the actual parameter is user-defined. It is recommended for ASVs calculated per sample (it may be ok to deactivate it for pooled sets)
df <- df[order(rowSums(df),decreasing=T),] # sort input by most abundant first (ignore first and last column)
# df <- df[,c(1:9,11:21,23:24)] # UPDATE 2020-02-28: This line was added to filter samples IL4 and HL4 and should be commented for the code to work correctly in any other set
if(is.na(group_sizes[1])){ # Update 2021-03-18: Added a default "no groups option where group size is set to nrow"
 	group_sizes=ncol(df)
 	groups=1
 	print(paste0("No groups provided, set default to 1 group with ", group_sizes," samples"))
}
cutoff1 <- (sum(df)*min_freq)
paste("The min freq of", min_freq, "corresponds to", cutoff1, "reads after removal of sample-singleton and items with <", min_reads, "reads")
df <- df[rowSums(df)>=(sum(df)*min_freq),] # Whole table filter features with less than min_freq total
bin_df <- df; bin_df[bin_df>0] = 1 #Copy and create a binary table (absense/presence); ignore first and last column)
init <- 0
compare <- data.frame(matrix(NA,nrow=nrow(df),ncol=groups)) # Create an empty matrix with the same number of otus with columns
if(groups==1){compare[,1]=rowSums(bin_df);compare=compare>=min_samples*group_sizes;}else{
	for(i in 1:groups){ # For each group, sum the columns in the binary table to compare them with the min_samples cutoff
		end <- init+group_sizes[i] #set the advancing window for columns
		compare[i] <- as.numeric(rowSums(bin_df[,(init+1):end])>=floor(group_sizes[i]*min_samples)) # Create a binary vector for results
		init <- init+group_sizes[i] #reset the starting column for next iteration
	}
}
df <- df[rowSums(compare)>=1,] # Filter the original set based on columns where in at least one group the feature passes the filter
df <- df[,colSums(df)>0] # Remove empty columns (if any)
print(groups)
write.table(df,output_name, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # and export it as a tsv file
