Actual commands: Started 2021-02-04
 #############################################
 ###        START ANALYSES HERE:           ###
 #############################################
# 1.- Input files
# We now have a total of 6 tables to analyze, we want to compare Mazatlan vs Lajitas ponds in an organ-based manner. To achieve this, we will first create sub graphs for each I and H organs:
cd /home/rod/Documents/02_Collaborations/Geneticas
ls 01_input_tables/genetics_lvl*.tsv
01_input_tables/genetics_lvl2.tsv  01_input_tables/genetics_lvl5.tsv  01_input_tables/genetics_lvlOTU.tsv
01_input_tables/genetics_lvl3.tsv  01_input_tables/genetics_lvl6.tsv
01_input_tables/genetics_lvl4.tsv  01_input_tables/genetics_lvl7.tsv

# 2.- Create split tables
mkdir 02_split_tables
# for organ in H I; do for lvl in 2 3 4 5 6 7 OTU;do echo 'cat 01_input_tables/genetics_lvl'$lvl'.tsv|Rscript Extract_groups_from_table.R 02_split_tables/genetics_lvl'$lvl' '$organ''; for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'_grp_'$organ'.tsv|Rscript Extract_groups_from_table.R 02_split_tables/genetics_lvl'$lvl'-'$organ' '$pond'';done;echo 'rm 02_split_tables/genetics_lvl'$lvl'_grp_'$organ'.tsv';done;done
for organ in H I; do for lvl in 2 3 4 5 6 7 OTU;do echo 'cat 01_input_tables/genetics_lvl'$lvl'.tsv|Rscript Extract_groups_from_table.R 02_split_tables/genetics_lvl'$lvl' '$organ''; for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'_grp_'$organ'.tsv|Rscript Extract_groups_from_table.R 02_split_tables/genetics_lvl'$lvl'-'$organ' '$pond'';done;done;done # same thing but avoid deleting whole H and I sets

# 3.- Calculate correlations
mkdir 03_corr_tables
for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Correlations.R 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond' pearson rows';done;done;done
for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do echo 'cat 02_split_tables/genetics_lvl'$lvl'_grp_'$organ'.tsv|Rscript Correlations.R 03_corr_tables/genetics_lvl'$lvl'-'$organ' pearson rows';done;done
# NOTE: We only use pearson correlations because spearman uses ranks and we don't care about ranks but rather exact numbers here

# 4.- Network analysis
mkdir 04_networks
cut -f 1 03_corr_tables/genetics_lvl*|sort|uniq|awk '{print $0 "\t" "A"}' >01_input_tables/meta_all.tsv # Create a dummy file containing all possible row names in all correlation tables we're working with. This was used as an alternative to creating one metadata file per table (there are 28 of them in total)

i=0.5; for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/meta_all.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows TRUE FALSE '$i'';done;done;done
# Also, for the whole talbes (H and I)
i=0.5; for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'-corr_pearson_rows.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/meta_all.tsv 04_networks/genetics_lvl'$lvl'-'$organ'-corr_pearson_rows TRUE FALSE '$i'';done;done


# 5.- Compare tables
mkdir 05_collated_comparisons

# The commands are in the script compare_centrality.R

Started 2021-02-11
 #############################################
 ###        Alternative A:           ###
 #############################################
# Instead of using the complete corr tables, we will filter them by the p value associated to correlation tests, then, we will use a 0.8 cutoff and use only OTU tables. Comparison procedures will remain unaltered but will be plotted with boxplot or violin plots.

mkdir 03_corr_tables
# We can safely resume from step 03 onward, using a different script that includes pvalue filtering
for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Correlations-filt_by_pval.R 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond' pearson rows 0.001';done;done;done|bash
# also with 0.01, because we kept way too few items with 0.001
# for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Correlations-filt_by_pval.R 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond' pearson rows 0.01';done;done;done|bash

# I got a three-column xref table linking OTUs with their respective total sum (relative) and taxonomy @ (/home/rod/Documents/02_Collaborations/Geneticas/01_input_tables/xref-OTUs_taxonomy.tsv)
# With this, we can create a simple phyla grouping:
cat 01_input_tables/xref-OTUs_taxonomy.tsv|sed -e 's/k__Bacteria; p__//' -e 's/; c__.*//' -e 's/taxonomy/Phylum/'|awk -F $'\t' 'BEGIN {OFS = FS}{print $1,$3,$2}' >01_input_tables/xref-OTUs_phyla.tsv
# For the taxa, we create an analogous table
for lvl in 2 3 4 5 6 7; do paste <(cut -f 1 01_input_tables/genetics_lvl$lvl.tsv) <(sed -e 's/\t.*//' -e 's/.*;p__//' -e 's/;.*//' 01_input_tables/genetics_lvl$lvl.tsv) <(awk -F'\t' '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' 01_input_tables/genetics_lvl$lvl.tsv);done| grep -v "#" >01_input_tables/xref-taxa_phyla.tsv
# And create a single table (in case it's required)
cat 01_input_tables/xref-OTUs_phyla.tsv 01_input_tables/xref-taxa_phyla.tsv >01_input_tables/xref-all-lvls_phyla.tsv
# The resulting metadata file has columns: #OTU ID Phylum  Rel_Sum

# Next, create the networks
mkdir 04_networks
i=0.8; for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows_maxpval_0.001.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/xref-all-lvls_phyla.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows-0.001 TRUE FALSE '$i'';done;done;done|bash

# Repeat with spearman
for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Correlations-filt_by_pval.R 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond' spearman rows 0.001';done;done;done|bash
i=0.8; for lvl in 2 3 4 5 6 7 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_spearman_rows_maxpval_0.001.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/xref-all-lvls_phyla.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_spearman_rows-0.001 TRUE FALSE '$i'';done;done;done|bash


# 5.- Compare tables
mkdir 05_collated_comparisons
# The remaining methods and plotting commandas are found in file /home/rod/Documents/02_Collaborations/Geneticas/compare_centrality.R 

# I modified the program to create an additional table depicting the composition of each subgraph larger than 10 nodes. Files are named with the -subgraphs.tsv prefix. By using this I can determine which taxa match each OTU since we have the actual OTU identifier and a cross-reference table of assigned taxonomy of each item (01_input_tables/xref-OTUs_taxonomy.tsv)

for organ in H I;do for grp in L A; do for i in $(cut -f 3 04_networks/genetics_lvlOTU-$organ\_grp_$grp\-corr_pearson_rows-0.001-cut-0.8-full-subgraphs.tsv);do paste <(grep -m 1 -w $i 04_networks/genetics_lvlOTU-$organ\_grp_$grp\-corr_pearson_rows-0.001-cut-0.8-full-subgraphs.tsv) <(grep -m 1 -w $i 01_input_tables/xref-OTUs_taxonomy.tsv| cut -f 2-);done >04_networks/subgraphs_lvlOTU-$organ\_grp_$grp\.tsv; done;done

# The resulting files were searched for probiotics from Pablo's work with no avail:
cut -f 2 01.1_input_files/probiotic80.txt|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_L.tsv 
cut -f 2 01.1_input_files/probiotic80.txt|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_A.tsv 
cut -f 2 01.1_input_files/probiotic80.txt|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_L.tsv 
cut -f 2 01.1_input_files/probiotic80.txt|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_A.tsv

# Repeat with the genus only (the set has few valid species-lvl identifications) but only vibrios, clostridia and shewanellas were found
cut -f 2 01.1_input_files/probiotic80.txt|sed -e 's/; s__.*//'|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_L.tsv # Only vibrios, clostridia
cut -f 2 01.1_input_files/probiotic80.txt|sed -e 's/; s__.*//'|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_A.tsv  # Only vibrios, clostridia and shewanellas
cut -f 2 01.1_input_files/probiotic80.txt|sed -e 's/; s__.*//'|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_L.tsv # Only vibrios, shewanellas
cut -f 2 01.1_input_files/probiotic80.txt|sed -e 's/; s__.*//'|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_A.tsv

# Now we try with Luigui's differential taxa (multiple levels)
sed 's/_/__/' 01.1_input_files/Diff_taxa_I_all_lvl-L.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_L.tsv >04_networks/Diff_taxa-I_L.tsv
sed 's/_/__/' 01.1_input_files/Diff_taxa_I_all_lvl-A.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_A.tsv >04_networks/Diff_taxa-I_A.tsv
sed 's/_/__/' 01.1_input_files/Diff_taxa_H_all_lvl-L.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_L.tsv >04_networks/Diff_taxa-H_L.tsv
sed 's/_/__/' 01.1_input_files/Diff_taxa_H_all_lvl-A.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_A.tsv >04_networks/Diff_taxa-H_A.tsv

sed 's/_/__/' 01.1_input_files/Diff_taxa_I_all_lvl-L.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_L.tsv
sed 's/_/__/' 01.1_input_files/Diff_taxa_I_all_lvl-A.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-I_grp_A.tsv
sed 's/_/__/' 01.1_input_files/Diff_taxa_H_all_lvl-L.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_L.tsv
sed 's/_/__/' 01.1_input_files/Diff_taxa_H_all_lvl-A.tsv|grep -wFf - 04_networks/subgraphs_lvlOTU-H_grp_A.tsv

 #############################################
 ###        Alternative B:           ###
 #############################################

# Now create the same type of table but only from phyla to genus (This cannot be created for the lower levels)
cat 01_input_tables/xref-OTUs_taxonomy.tsv|sed -e 's/k__Bacteria; //' -e 's/; s__.*//'|awk -F $'\t' 'BEGIN {OFS = FS}{print $1,$3,$2}' >01_input_tables/xref-OTUs_phlyla-genera.tsv
# And one for the lvl 6 (genus) tables
for lvl in 6; do paste <(cut -f 1 01_input_tables/genetics_lvl$lvl.tsv) <(sed -e 's/\t.*//' -e 's/k__Bacteria;//' -e 's/s__.*//' 01_input_tables/genetics_lvl$lvl.tsv) <(awk -F'\t' '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' 01_input_tables/genetics_lvl$lvl.tsv);done| grep -v "#" >01_input_tables/xref-lvl6_phyla-genera.tsv
cat 01_input_tables/xref-OTUs_phlyla-genera.tsv 01_input_tables/xref-lvl6_phyla-genera.tsv >01_input_tables/xref_lvl6nOTU_phyla-genera.tsv

mkdir 04_networks
i=0.8; for lvl in 6 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows_maxpval_0.001.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/xref_lvl6nOTU_phyla-genera.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows-0.001 TRUE FALSE '$i'';done;done;done

# There were way too many items, we'll use only family instead

 #############################################
 ###        Alternative C:           ###
 #############################################

# Now create the same type of table but only from family only (This cannot be created for the lower levels)
cat 01_input_tables/xref-OTUs_taxonomy.tsv|sed -e 's/k__.*f__//' -e 's/;.*//' -e 's/\t$/\tUndefined/'|awk -F $'\t' 'BEGIN {OFS = FS}{print $1,$3,$2}' >01_input_tables/xref-OTUs_families.tsv

# And one for the lvl 6 (genus) tables containing families as well
for lvl in 6; do paste <(cut -f 1 01_input_tables/genetics_lvl$lvl.tsv) <(sed -e 's/\t.*//' -e 's/k__.*f__//' -e 's/;.*//' -e 's/^$/Undefined/' 01_input_tables/genetics_lvl$lvl.tsv) <(awk -F'\t' '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' 01_input_tables/genetics_lvl$lvl.tsv);done| grep -v "#" >01_input_tables/xref-lvl6_families.tsv
cat 01_input_tables/xref-OTUs_families.tsv 01_input_tables/xref-lvl6_families.tsv >01_input_tables/xref_lvl6nOTU_families.tsv

mkdir 04_networks
i=0.8; for lvl in 6 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows_maxpval_0.001.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/xref_lvl6nOTU_families.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows-0.001 TRUE FALSE '$i'';done;done;done

 #############################################
 ###        Alternative D:           ###
 #############################################
# Instead of using the complete corr tables, we will filter them by the p value associated to correlation tests, then, we will use a 0.7 cutoff and pval 0.01 and use only OTU and lvl 6 tables. Comparison procedures will remain unaltered but will be plotted with boxplot or violin plots.

mkdir 03_corr_tables
# We can safely resume from step 03 onward, using a different script that includes pvalue filtering
# also with 0.01, because we kept way too few items with 0.001
for lvl in 6 OTU; do for organ in H I; do for pond in L A; do echo 'cat 02_split_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Correlations-filt_by_pval.R 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond' pearson rows 0.01';done;done;done|bash

# Next, create the networks
mkdir 04_networks
i=0.7; for lvl in 6 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows_maxpval_0.01.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/xref-all-lvls_phyla.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows-0.01 TRUE FALSE '$i'';done;done;done|bash

# 5.- Compare tables
mkdir 05_collated_comparisons
# The remaining methods and plotting commandas are found in file /home/rod/Documents/02_Collaborations/Geneticas/compare_centrality.R 


 #############################################
 ###        Alternative E:           ###
 #############################################
# We will now test what happens when we consider only items with over 50% items.
# IMPORTANT NOTE: There is a problem with the cutoff, since group IA has 7 items. We have decided to remove sample IA9, column 24  (the last one but this is still arbitrary)

# 2.- Create split tables
mkdir 02_split_tables-noIA9
for organ in H I; do for lvl in 2 3 4 5 6 7 OTU;do echo 'cat 01_input_tables/genetics_lvl'$lvl'.tsv|cut -f 1-23,25-|Rscript Extract_groups_from_table.R 02_split_tables-noIA9/genetics_lvl'$lvl' '$organ''; for pond in L A; do echo 'cat 02_split_tables-noIA9/genetics_lvl'$lvl'_grp_'$organ'.tsv|Rscript Extract_groups_from_table.R 02_split_tables-noIA9/genetics_lvl'$lvl'-'$organ' '$pond'';done;done;done|bash

# We use the same split tables in 02_split_tables and script Sparsity_reduction_merged_tables.R
mkdir 02.1_split_tables-50perc_of_samples
for organ in H I; do for lvl in 2 3 4 5 6 7 OTU;do for pond in L A; do echo 'cat 02_split_tables-noIA9/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Sparsity_reduction_merged_tables.R 0 0 0.5 02.1_split_tables-50perc_of_samples/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv';done;done;done

mkdir 03_corr_tables
# We can safely resume from step 03 onward, using a different script that includes pvalue filtering
# also with 0.01, because we kept way too few items with 0.001
for lvl in 6 OTU; do for organ in H I; do for pond in L A; do echo 'cat 02_split_tables-noIA9/genetics_lvl'$lvl'-'$organ'_grp_'$pond'.tsv|Rscript Correlations-filt_by_pval.R 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond' pearson rows 0.01';done;done;done|bash

# Next, create the networks
mkdir 04_networks
i=0.7; for lvl in 6 OTU; do for organ in H I; do for pond in L A; do echo 'cat 03_corr_tables/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows_maxpval_0.01.tsv|Rscript Network_from_symmetrical_matrix.R 01_input_tables/xref-all-lvls_phyla.tsv 04_networks/genetics_lvl'$lvl'-'$organ'_grp_'$pond'-corr_pearson_rows-0.01 TRUE FALSE '$i'';done;done;done|bash

# 5.- Compare tables
mkdir 05_collated_comparisons
# The remaining methods and plotting commandas are found in file /home/rod/Documents/02_Collaborations/Geneticas/compare_centrality.R 
