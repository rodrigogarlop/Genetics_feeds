# This will be a fixed script that will only work for plotting the density plot of the distributions comparing A and L ponds for each organ.
# We want to compare all centrality measurements at each tax lvl
# It won't be usefull for any other comparison.
setwd("/home/rod/Documents/02_Collaborations/Geneticas")
lvl <- c("2","3","4","5","6","7","OTU")
# tab1 <- "04_networks/genetics_lvl5-H_grp_A-corr_pearson_rows-cut-0.5-full-node_specs.tsv"
# tab2 <- "04_networks/genetics_lvl5-H_grp_L-corr_pearson_rows-cut-0.5-full-node_specs.tsv"
compare_tables <- function(tab1, tab2, column){ # Use this to go through each pair of files
	df1 <- read.table(tab1, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # First, load each file
	df2 <- read.table(tab2, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
	df1 <- df1[1:(nrow(df1)-2),] # and remove global features
	df2 <- df2[1:(nrow(df2)-2),]
	rownames(df1) <- df1[,1] # also, vertices names (v1...n don't mean anything), so replace them with actual taxa/cluster names
	rownames(df2) <- df2[,1]
	temp <- c(rownames(df1),rownames(df2))
	temp <- temp[duplicated(temp)] # next, keep only matching taxa/clusters and in the same order (well use the first)
	df1 <- df1[temp,] # Now we have two alphabetically sorted, matching-item tables from both files
	df2 <- df2[temp,]
	out <- cbind(df1[column],df2[column]) # Extract the target column as data.frame
	return(out)
}
transp_cols <- function(color, percent = 50, name = NULL) { # Set transparent color vector: Function from ## www.dataanalytics.org.uk
	rgb.val <- col2rgb(color)
	t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100 - percent) * 255 / 100, names = name)
	invisible(t.col) ## Save the color
}
plot_density <- function(tab, lvl){ # from a 2-column table, plot densities
	a <- tab[,1][!is.na(tab[,1])]
	b <- tab[,2][!is.na(tab[,2])]
	A <- stats::density(a)
	B <- stats::density(b)
	ab <- ((!is.na(tab[,1])))*((!is.na(tab[,2]))) # determine which are directly comparable (ocurring in 2 gprs)
	range.x <- range(c(A[[1]],B[[1]]))
	max.y <- max(c(A[[2]],B[[2]]))
	pval <- wilcox.test(tab[as.logical(ab),1],tab[as.logical(ab),2], paired=T)$p.value
	pval2 <- t.test(tab[as.logical(ab),1],tab[as.logical(ab),2], paired=T)$p.value
	if(pval < 0.0001){pval = "< 0.0001"}else{pval=round(pval,4)}
	if(pval2 < 0.0001){pval2 = "< 0.0001"}else{pval2=round(pval2,4)}
	cor <- round(cor(tab[as.logical(ab),1],tab[as.logical(ab),2]),4)
	plot(A, xlim=range.x, ylim=c(0,max.y), yaxt='n', xlab="", ylab=paste0("Lvl ",lvl), main="") # Create the density graph
	polygon(A,col=t_cols[1],border=solid_cols[1]) # and color it accordingly
	polygon(B,col=t_cols[2],border=solid_cols[2]) # and color it accordingly
	legend("topleft", legend=(c(paste0("Wilcoxon p-value: ",pval), paste0("t p-value: ",pval2), paste0("Pearson corr: ", cor))), bty = "n")
	plot(A,B)
}
process_column <- function(column=2,organ){ # this is just to circle through columns and plot
	for(i in lvl){
		tab1 <- paste0("04_networks/genetics_lvl",i,"-",organ,"_grp_A-corr_pearson_rows-cut-0.5-full-node_specs.tsv")
		tab2 <- paste0("04_networks/genetics_lvl",i,"-",organ,"_grp_L-corr_pearson_rows-cut-0.5-full-node_specs.tsv")
		tab <- compare_tables(tab1, tab2, column)
		plot_density(tab,i)
# 		plot(tab)
	}
}

 ### Main ###
t_cols <- c(transp_cols("cornflowerblue",50),transp_cols("firebrick",50))
solid_cols <- c("blue","red")
pdf("05_collated_comparisons/degree_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(2,"H"); dev.off()
pdf("05_collated_comparisons/degree_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(2,"I"); dev.off()
pdf("05_collated_comparisons/eigenvector_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(3,"H"); dev.off()
pdf("05_collated_comparisons/eigenvector_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(3,"I"); dev.off()
pdf("05_collated_comparisons/closeness_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(4,"H"); dev.off()
pdf("05_collated_comparisons/closeness_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(4,"I"); dev.off()
pdf("05_collated_comparisons/betweenness_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(5,"H"); dev.off()
pdf("05_collated_comparisons/betweenness_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(5,"I"); dev.off()
pdf("05_collated_comparisons/transitivity_u-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(6,"H"); dev.off()
pdf("05_collated_comparisons/transitivity_u-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(6,"I"); dev.off()
pdf("05_collated_comparisons/closeness_w-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(7,"H"); dev.off()
pdf("05_collated_comparisons/closeness_w-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(7,"I"); dev.off()
pdf("05_collated_comparisons/betweenness_w-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(8,"H"); dev.off()
pdf("05_collated_comparisons/betweenness_w-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(8,"I"); dev.off()
pdf("05_collated_comparisons/hub_score-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(9,"H"); dev.off()
pdf("05_collated_comparisons/hub_score-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(9,"I"); dev.off()
pdf("05_collated_comparisons/vertex_strength-H.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(10,"H"); dev.off()
pdf("05_collated_comparisons/vertex_strength-I.pdf"); par(mfrow = c(7,1), mar = c(2, 4, 0.2, 0.5)); process_column(10,"I"); dev.off()




# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)
# plot_density(tab)

# Alternative A (this was ultimately used), we'll analyze only the OTU sets
compare_tables <- function(tab1, tab2, column){ # Use this to go through each pair of files
	df1 <- read.table(tab1, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # First, load each file
	df2 <- read.table(tab2, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
	df1 <- df1[1:(nrow(df1)-2),] # and remove global features
	df2 <- df2[1:(nrow(df2)-2),]
	rownames(df1) <- df1[,1] # also, vertices names (v1...n don't mean anything), so replace them with actual taxa/cluster names
	rownames(df2) <- df2[,1]
	temp <- c(rownames(df1),rownames(df2))
	temp <- temp[duplicated(temp)] # next, keep only matching taxa/clusters and in the same order (well use the first)
	df1 <- df1[temp,] # Now we have two alphabetically sorted, matching-item tables from both files
	df2 <- df2[temp,]
	out <- cbind(df1[column],df2[column]) # Extract the target column as data.frame
	return(out)
}
transp_cols <- function(color, percent = 50, name = NULL) { # Set transparent color vector: Function from ## www.dataanalytics.org.uk
	rgb.val <- col2rgb(color)
	t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100 - percent) * 255 / 100, names = name)
	invisible(t.col) ## Save the color
}
plot_density <- function(tab, lvl){ # from a 2-column table, plot densities
	a <- tab[,1][!is.na(tab[,1])]
	b <- tab[,2][!is.na(tab[,2])]
	A <- stats::density(a)
	B <- stats::density(b)
	ab <- ((!is.na(tab[,1])))*((!is.na(tab[,2]))) # determine which are directly comparable (ocurring in 2 gprs)
	range.x <- range(c(A[[1]],B[[1]]))
	max.y <- max(c(A[[2]],B[[2]]))
	pval <- wilcox.test(tab[as.logical(ab),1],tab[as.logical(ab),2], paired=T)$p.value
	pval2 <- t.test(tab[as.logical(ab),1],tab[as.logical(ab),2], paired=T)$p.value
	if(pval < 0.0001){pval = "< 0.0001"}else{pval=round(pval,4)}
	if(pval2 < 0.0001){pval2 = "< 0.0001"}else{pval2=round(pval2,4)}
	cor <- round(cor(tab[as.logical(ab),1],tab[as.logical(ab),2]),4)
	plot(A, xlim=range.x, ylim=c(0,max.y), yaxt='n', xlab="", ylab=paste0("Lvl ",lvl), main=colnames(tab)[1]) # Create the density graph
	polygon(A,col=t_cols[1],border=solid_cols[1]) # and color it accordingly
	polygon(B,col=t_cols[2],border=solid_cols[2]) # and color it accordingly
	legend("topleft", legend=(c(paste0("Wilcoxon p-value: ",pval), paste0("t p-value: ",pval2), paste0("Pearson corr: ", cor))), bty = "n")
# 	plot(A,B)
}
plot_boxplot <- function(tab, lvl){ # from a 2-column table, plot densities
	a <- tab[,1][!is.na(tab[,1])]
	b <- tab[,2][!is.na(tab[,2])]
# 	range.x <- range(c(A[[1]],B[[1]]))
# 	max.y <- max(c(A[[2]],B[[2]]))
	pval <- wilcox.test(a,b, paired=F)$p.value
	pval2 <- t.test(a,b, paired=F)$p.value
	if(pval < 0.0001){pval = "< 0.0001"}else{pval=round(pval,4)}
	if(pval2 < 0.0001){pval2 = "< 0.0001"}else{pval2=round(pval2,4)}
	boxplot(a,b, col=t_cols, frame.plot=F, border=solid_cols, main=names(tab)[1], xaxt='n', outline=T)
	axis(1,las=1, labels=c("Lajitas","Mazatlan"), at=1:2)
	legend("topright", legend=(c(paste0("Mann-Whitney p-value: ",pval), paste0("t p-value: ",pval2))), bty = "n")
# 	plot(A,B)
}
process_column <- function(column=2,organ){ # this is just to circle through columns and plot
	for(i in lvl){
		tab1 <- paste0("04_networks/genetics_lvl",i,"-",organ,"_grp_A-corr_pearson_rows-0.01-cut-0.7-full-node_specs.tsv")
		tab2 <- paste0("04_networks/genetics_lvl",i,"-",organ,"_grp_L-corr_pearson_rows-0.01-cut-0.7-full-node_specs.tsv")
		tab <- compare_tables(tab1, tab2, column)
		plot_density(tab,i)
		plot_boxplot(tab)
	}
}
lvl="OTU"
t_cols <- c(transp_cols("cornflowerblue",50),transp_cols("firebrick",50))
solid_cols <- c("blue","red")
pdf(paste0("05_collated_comparisons/lvl",lvl,"-degree_u-H.pdf")); par(las=2); process_column(2,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-degree_u-I.pdf")); par(las=2); process_column(2,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-eigenvector_u-H.pdf")); par(las=2); process_column(3,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-eigenvector_u-I.pdf")); par(las=2); process_column(3,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-closeness_u-H.pdf")); par(las=2); process_column(4,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-closeness_u-I.pdf")); par(las=2); process_column(4,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-betweenness_u-H.pdf")); par(las=2); process_column(5,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-betweenness_u-I.pdf")); par(las=2); process_column(5,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-transitivity_u-H.pdf")); par(las=2); process_column(6,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-transitivity_u-I.pdf")); par(las=2); process_column(6,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-closeness_w-H.pdf")); par(las=2); process_column(7,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-closeness_w-I.pdf")); par(las=2); process_column(7,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-betweenness_w-H.pdf")); par(las=2); process_column(8,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-betweenness_w-I.pdf")); par(las=2); process_column(8,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-hub.score-H.pdf")); par(las=2); process_column(9,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-hub.score-I.pdf")); par(las=2); process_column(9,"I"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-vertex.strength-H.pdf")); par(las=2); process_column(10,"H"); dev.off()
pdf(paste0("05_collated_comparisons/lvl",lvl,"-vertex.strength-I.pdf")); par(las=2); process_column(10,"I"); dev.off()

# [1] "vertex_name"     "degree_u"        "eigenvector_u"   "closeness_u"    
#  [5] "betweenness_u"   "transitivity_u"  "closeness_w"     "betweenness_w"  
#  [9] "hub.score"       "vertex.strength"
