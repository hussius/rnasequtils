require(calibrate)

# Plots like "Top X% have Y% of the tags"
# Store a vector for each gene

real.estate <- function(counts, plot=FALSE){
require(Hmisc)
tagFperProc = matrix(nrow = ncol(counts), ncol = 100)

for (j in 1:ncol(counts)){
x <- counts[order(counts[,j], decreasing=TRUE),j]
#print(length(which(x==0))) 

# Plots like "Top X% have Y% of the tags"
gFrac <- vector()
tFrac <- vector()
for (i in seq(1:100)){
	nGenes <- floor((i/100) * nrow(counts))
	nTags <- sum(x[1:nGenes])
	tTags <- sum(x)
	gFrac <- c(gFrac, i)
	tFrac <- c(tFrac, nTags/tTags)
}

tagFperProc[j,] <- tFrac

if(plot==TRUE){
	plot(tFrac, gFrac, xlab = "Fraction of total #sequences", ylab="Top % of ENSEMBL transcripts", xlim=c(0,1), main=tissue[j], ylim=c(0,20))
	lines(tFrac, gFrac)
	readline()
	}

}

meanF <- colMeans(tagFperProc)
stdF <- apply(tagFperProc, 2, sd)

meanPerc <- 100 * meanF
stdPerc <- 100 * stdF

if(plot==TRUE){
errbar(1:100, meanPerc, yplus = meanPerc + stdPerc, yminus = meanPerc - stdPerc, xlim=c(0,20), ylim=c(0,25), xlab = "Top X% of transcripts", ylab = "% of total sequences in top X% transcripts")

boxplot(tagFperProc, xlim=c(0,20), xlab = "Top X% of transcripts", ylab = "% of total sequences in top X% transcripts")
}

return(tagFperProc)
	
}
load.fpkm <- function(fname){
	fpkm <- read.delim(fname,stringsAsFactors=FALSE)
	return(fpkm)
}

load.counts <- function(fname){
	counts <- read.delim(fname,stringsAsFactors=FALSE)
	return(counts)
}

do.SVD = function(m, comp.1=1, comp.2=2){ # returns eig.cell
	s <- svd(m)
	ev <- s$d^2 / sum(s$d^2)
	return(s$u[,c(comp.1, comp.2)])
}

project.SVD <- function(m, eig.cell){
	return(t(m) %*% eig.cell)
}

plot.SVD <- function(m, comp.1=1, comp.2=2, groups=rep("blue", ncol(m))){
	eig <- do.SVD(m, comp.1, comp.2)
	proj <- project.SVD(m, eig)
	plot(proj, pch=20, col="white")
	points(proj, col=as.numeric(as.factor(groups)), pch=20)
	textxy(proj[,1],proj[,2],labs=colnames(m))
	return(proj)
}

loadings.SVD <- function(m, comp=1, gene.ids = rownames(m)){
	s <- svd(m)
	l <- s$u[,comp]
	names(l) <- gene.ids
	l.s <- l[order(l)]
	return(l.s)
}

plot.loadings.SVD <- function(m, comp=1, cutoff=0.1, gene.ids = rownames(m)){
	l <- loadings.SVD(m, comp, gene.ids)
	barplot(l[abs(l)>cutoff],las=2,main=paste("PC", comp, "cutoff", cutoff),cex.names=0.6)
}

compare.pearson <- function(matrix.list, labels=NULL){
	par(mfrow=c(1,length(matrix.list)))
	for (i in 1:length(matrix.list)){
		hist(cor(matrix.list[[i]]),n=100,main=labels[i])
	}
}

compare.spearman <- function(matrix.list,labels=NULL){
	par(mfrow=c(1,length(matrix.list)))
	for (i in 1:length(matrix.list)){
		hist(cor(matrix.list[[i]],method="spearman"),n=100,main=labels[i])
	}
}

arcsintrans.gene <- function(count.vec){
	return(sqrt(sum(count.vec))*asin(sqrt(count.vec/sum(count.vec))))	
}

arcsintrans <- function(m){
	return(apply(m,2,arcsintrans.gene))
}

high.outliers.vec <- function(vec, nosd = 4){
	return(vec > mean(vec)+nosd*sd(vec))
}

ntrue <- function(boolvec){
	return(length(which(boolvec==TRUE)))
}

count.outliers.per.sample <- function(m, no.sd = 4){
	o <- apply(m, 1, high.outliers.vec, nosd = no.sd)
	return(apply(o,1,ntrue))
}

get.outliers.for.sample <- function(m, samp, no.sd = 4){
	o <- apply(m, 1, high.outliers.vec, nosd = no.sd)
	olgenes <- names(which(o[samp,]==TRUE))
	return(olgenes)
}

plot.expr.from.list <- function(m, de, gene.names=NULL){
	if (is.null(gene.names)){
		for (g in de[,1]){
			barplot(as.matrix(m[g,],las=2,main=g))
			readline()
		}}
	else{
		for (g in de[,1]){
			barplot(as.matrix(m[gene.names==g,],las=2,main=g))
			readline()
		}
	}
}

normalize.voom <- function(counts){
	require(limma)
	return(voom(counts)$E)
}
 
cpm.tmm <- function(counts, groups=NA){
	require(edgeR)
	if(is.na(groups)){
		d<-DGEList(counts=counts)
	}
	else{
		d<-DGEList(counts=counts, group=groups)
	}
	d <- calcNormFactors(d, method="TMM") 
	return(cpm(d, normalized.lib.sizes=TRUE))
}

# Example usage

#counts <- load.counts("/Users/mikaelhuss/Desktop/HiSeq/t_olsson_11_02/count_table.txt")

#plot.SVD(counts)
#plot.SVD(counts, comp.1=2, comp.2=3)
#plot.loadings.SVD(counts)

# Need to treat the FPKM table differently
#fpkm <- load.fpkm("/Users/mikaelhuss/Desktop/HiSeq/t_olsson_11_02/fpkm_table.txt")
#f.num <- as.matrix(fpkm[,3:ncol(fpkm)])
#f.ens.id <- fpkm[,1]
#f.gene.id <- fpkm[,2]

# Now it should work to do everything on FPKMs
#plot.SVD(f.num)
#plot.loadings.SVD(f.num, gene.ids=f.gene.id)

# Arcsine transformation
#arc <- arcsintrans(counts)
#log2plus1 <- log(counts+1, base=2)
#compare.pearson(list(f.num, arc))
#compare.pearson(list(f.num, counts, arc, log2plus1),labels=c("FPKM","Raw counts","Arcsin","log2 (1+counts)"))

# Finding the samples with the most outliers, and the specific values that are outliers for a sample
#ol.counts <- count.outliers.per.sample(arc)

# Try to find the genes which are most commonly outliers (but does that make sense?)
# Maybe export this list and see if there is some sort of overrepresentation
#ol.genes <- vector()
#idx <- 1
#for(name in colnames(arc)){
#	ol.genes <- c(ol.genes, get.outliers.for.sample(arc, name))
#	print(idx)
#	idx <- 1 + idx
#}

#table(ol.genes) # <- which genes are most commonly outliers
#ol.u <- unique(ol.genes) # which unique genes have been outliers

# Plot differentially expressed genes
#de <- read.delim("/Users/mikaelhuss/Desktop/HiSeq/t_olsson_11_02/DE/bayseq-genes.txt", sep="\n", stringsAsFactors=F, header=T)

#plot.expr.from.list(f.num, de, gene.names=f.ens.id)
#plot.expr.from.list(arc, de)