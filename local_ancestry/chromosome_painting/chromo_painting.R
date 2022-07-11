# chromosome_painting_genomewide.R
# This script first identifies Ancestry-Informative MArkers (AIMs) in parental pools given some allele frequency difference threshold (hard-coded, thresh = .8), and then genotype AIMs in hybrids.

### libs, working directory, and stuff
options(stringsAsFactors = F)
library(stringr)

### Set which scaffold is analysed now
my.scaffold = commandArgs(trailingOnly=TRUE)[1]
setwd(paste0("/scratch/project_2001099/nouhaudp/reseq/snp/analyses/local_ancestry/chromopaint/", my.scaffold))


###
#### Read data ----
###

res.dir = paste0("/scratch/project_2001099/nouhaudp/reseq/snp/analyses/local_ancestry/chromopaint/", my.scaffold, "/")
geno = read.table("/scratch/project_2001099/nouhaudp/reseq/snp/geno/nonphased_out.geno.gz", h = T, comment.char = "")
names(geno)[1] = "CHROM"
geno = geno[geno$CHROM == my.scaffold,]
geno = geno[,-which(names(geno) == "OUT")] # removes the outgroup from the geno file
smpl_table = read.table("/scratch/project_2001099/nouhaudp/reseq/sample_table_R_ordered_LanFix.tab", h=T, comment.char = "")
smpl_table = smpl_table[smpl_table$sex == "f",]
fai = read.table("/scratch/project_2001099/nouhaudp/reseq/ref/Formica_hybrid_v1.fa.fai", h = FALSE)




###
#### Prep data ----
###

### 1. Parental diploid genotypes are haploidized to help with allele counts

print(paste0("### Extracting alleles for parental individuals..."))

# Order columns per species
geno.ord = geno[, c(1, 2, match(smpl_table$id, names(geno[-c(1, 2)]))+2)]

# Split diploid genotypes per species - aq outside Finland
Aq.idx = which(smpl_table$pop.merged == "Aq")+2
aq = geno.ord[, Aq.idx]
females1 = females2 = matrix(data = NA, ncol = ncol(aq), nrow = nrow(aq))
for (i in 1:ncol(aq)){
  females1[, i] = str_split_fixed(string = aq[,i], pattern = "/", n = 2)[,1]
  females2[, i] = str_split_fixed(string = aq[,i], pattern = "/", n = 2)[,2]
}
aq = cbind(females1, females2)[,rep(seq(1, ncol(aq)), each = 2)+c(0, ncol(aq))]
colnames(aq) = paste0(rep(smpl_table$id[Aq.idx-2], each = 2), "-", c(1, 2))

# Split diploid genotypes per species - pol outside Finland
Pol.idx = which(smpl_table$pop.merged == "Pol")+2
pol = geno.ord[, Pol.idx]
females1 = females2 = matrix(data = NA, ncol = ncol(pol), nrow = nrow(pol))
for (i in 1:ncol(pol)){
  females1[, i] = str_split_fixed(string = pol[,i], pattern = "/", n = 2)[,1]
  females2[, i] = str_split_fixed(string = pol[,i], pattern = "/", n = 2)[,2]
}
pol = cbind(females1, females2)[,rep(seq(1, ncol(pol)), each = 2)+c(0, ncol(pol))]
colnames(pol) = paste0(rep(smpl_table$id[Pol.idx-2], each = 2), "-", c(1, 2))

# combine both species & tag missing data
all.parent = cbind(aq, pol)
all.parent[all.parent == "N"] = NA


### 2. Filtering ----

print(paste0("### Filtering..."))

# Flag sites with more than 12.5% missing data overall in both referecne panels (3/((6*2)*2))
NA.sites = which(rowSums(is.na(all.parent)) > 3)

# Flag all non-biallelic sites (incl. triallelic, but none at this step in theory)
monomorphic.sites = apply(all.parent, 1, function(x) length(na.omit(unique(x))))
monomorphic.sites = which(monomorphic.sites != 2)

# Remove sites
removed.sites = unique(sort(c(NA.sites, monomorphic.sites)))
print(paste0("# Number of sites with more than 12.5% missing data: ", length(NA.sites)))
print(paste0("# Number of monomorphic sites: ", length(monomorphic.sites)))
print(paste0("# Total number of sites removed (fraction): ", length(removed.sites),
             " (", round(length(removed.sites)/nrow(all.parent)*100, 2), "%)"))
all.parent = all.parent[-removed.sites,]
print(paste0("# Remaining number of sites: ", nrow(all.parent)))


#### 3. Count parental alleles ----

# count alleles
print(paste0("### Counting parental alleles..."))
parent.cnts = data.frame(CHR = geno.ord$CHROM[-removed.sites],
                         POS = geno.ord$POS[-removed.sites],
                         allele1 = rep(NA, nrow(all.parent)),
                         allele2 = rep(NA, nrow(all.parent)),
                         cnts.a1.aq = rep(NA, nrow(all.parent)),
                         cnts.a2.aq = rep(NA, nrow(all.parent)),
                         cnts.N.aq = rep(NA, nrow(all.parent)),
                         cnts.a1.pol = rep(NA, nrow(all.parent)),
                         cnts.a2.pol = rep(NA, nrow(all.parent)),
                         cnts.N.pol = rep(NA, nrow(all.parent)))

for (i in 1:nrow(all.parent)){
  
  print(paste(round(i/nrow(all.parent), 3)*100, "%"))
  
  # alleles
  my.alleles = names(sort(table(all.parent[i,]), decreasing = TRUE))
  parent.cnts[i, 3:4] = my.alleles
  
  # aq counts
  parent.cnts[i, 5] = sum(all.parent[i,1:ncol(aq)] == parent.cnts[i, 3], na.rm = TRUE)
  parent.cnts[i, 6] = sum(all.parent[i,1:ncol(aq)] == parent.cnts[i, 4], na.rm = TRUE)
  
  # pol counts
  parent.cnts[i, 8] = sum(all.parent[i,(1:ncol(pol))+ncol(aq)] == parent.cnts[i, 3], na.rm = TRUE)
  parent.cnts[i, 9] = sum(all.parent[i,(1:ncol(pol))+ncol(aq)] == parent.cnts[i, 4], na.rm = TRUE)
  
}

# missing sites
parent.cnts$cnts.N.aq = ncol(aq) - rowSums(parent.cnts[,5:6])
parent.cnts$cnts.N.pol = ncol(pol) - rowSums(parent.cnts[,8:9])

# compute a1 frequency in each group
parent.cnts$freq.a1.aq = parent.cnts$cnts.a1.aq/rowSums(parent.cnts[,5:6])
parent.cnts$freq.a1.pol = parent.cnts$cnts.a1.pol/rowSums(parent.cnts[,8:9])

# compute absolute frequency difference with a1 counts
parent.cnts$abs.diff.a1 = abs(parent.cnts$freq.a1.aq - parent.cnts$freq.a1.pol)


# Save file
write.table(parent.cnts,
            paste0(res.dir, "allele_counts_all_parents_noFin_", my.scaffold, "_biall_125percNA.tsv"),
            col.names = T, row.names = F, quote = F, sep = "\t")

# If needed / debugging
# parent.cnts = read.table(paste0(res.dir, "allele_counts_all_parents_noFin_allScafs_biall_125percNA.tsv"), h=T)



#### 4. Identify ancestry-informative markers AIMs ----

# cumulative numbers
n.haploid.samples = rowSums(parent.cnts[1, 5:7])
freq.diff.summary = data.frame(cutoff = seq(0, 1, by = 1/n.haploid.samples),
                               n = rep(NA, n.haploid.samples+1))

for (i in 1:nrow(freq.diff.summary)){
  freq.diff.summary$n[i] = sum(round(parent.cnts$abs.diff.a1, 2) == round(freq.diff.summary$cutoff[i],2))
}

freq.diff.summary$revcum.n = rev(cumsum(rev(freq.diff.summary$n)))

write.table(freq.diff.summary,
  paste0(res.dir, "AFD_all_parents_noFin_", my.scaffold, "_biall_125percNA.tsv"),
  col.names = T, row.names = F, quote = F, sep = "\t")




#### 5. Extract AIMs ----

thresh = .8
print(paste0("### Identifying AIMs with frequencies differing from at least ", thresh*100, "%..."))

pos.aim = paste0(parent.cnts$CHR[which(round(parent.cnts$abs.diff.a1, 2) >= thresh)],
                 "_",
                 parent.cnts$POS[which(round(parent.cnts$abs.diff.a1, 2) >= thresh)])
aq.alleles.aim = parent.cnts[which(round(parent.cnts$abs.diff.a1, 2) >= thresh),]
summary(aq.alleles.aim)

cov = aq.alleles.aim[,5:6]
cov.max = apply(cov, 1, which.max)
al = aq.alleles.aim[,3:4]

aq.alleles.aim = NULL

for(i in 1:nrow(al)){
  aq.alleles.aim = c(aq.alleles.aim, al[i, cov.max[i]])
}

head(aq.alleles.aim)


#### 6. Prep genotypes for hybrids & parents ----

print(paste0("### Prepare hybrid genotypes for painting..."))

# Subset AIMs
tmp.pos = paste0(geno.ord$CHROM, "_", geno.ord$POS)
dat = geno.ord[which(tmp.pos %in% pos.aim), ]

# Split genotypes
hyb = dat[, -c(1:2)]
hyb.chr1 = hyb.chr2 = matrix(data = NA, ncol = ncol(hyb), nrow = nrow(hyb))
for (i in 1:ncol(hyb)){
    hyb.chr1[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,1]
    hyb.chr2[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,2]
  }
}

colnames(hyb.chr1) = colnames(hyb.chr2) = names(hyb)

# Tag missing data
hyb.chr1[hyb.chr1 == "N"] = NA
hyb.chr2[hyb.chr2 == "N"] = NA

# Male filtering | deprecated
hyb.chr1.noNAmales = hyb.chr1 #[-which(rowSums(is.na(hyb.chr1))==40),]
hyb.chr2.noNAmales = hyb.chr2 #[-which(rowSums(is.na(hyb.chr1))==40),] # want to remove the same sites
aq.alleles.aim.noNAmales = aq.alleles.aim #[-which(rowSums(is.na(hyb.chr1))==40)]
pos.aim.noNAmales = pos.aim #[-which(rowSums(is.na(hyb.chr1))==40)]



#### 7. Compute HI for hybrids & parents ----

print(paste0("### Computing HIs..."))

# Sum over sites
head(aq.alleles.aim.noNAmales)
n.aq.alleles.chr1 = colSums(hyb.chr1.noNAmales == aq.alleles.aim.noNAmales, na.rm = T)
n.aq.alleles.chr2 = colSums(hyb.chr2.noNAmales == aq.alleles.aim.noNAmales, na.rm = T)
n.tot.alleles = (2*nrow(hyb.chr1.noNAmales))

# Compute HI
HI.aim = (n.aq.alleles.chr1 + n.aq.alleles.chr2)/n.tot.alleles
print(HI.aim)



#### Chromosome painting ----

print(paste0("### Painting..."))

real.pos.aim.noNAmales = as.numeric(matrix(ncol = 2, byrow = T,
                                           data = unlist(strsplit(x = pos.aim.noNAmales, split = "_")))[, 2])
tmp1 = hyb.chr1.noNAmales == aq.alleles.aim.noNAmales
tmp2 = hyb.chr2.noNAmales == aq.alleles.aim.noNAmales
paint.dat = tmp1+tmp2

# save painting data
writedat = cbind(rep(my.scaffold, length(real.pos.aim.noNAmales)), real.pos.aim.noNAmales, paint.dat)
colnames(writedat)[1:2] = c("CHR", "POS")
write.table(x = writedat,
            file = paste0("painted_genotypes_all_", thresh*100, "percAFD_", my.scaffold, ".tsv"),
            col.names = T, row.names = F, quote = F, sep = "\t")

# Create a vector of y positions for each haplotype
my.y = smpl_table$pop.idx+c(0:59)
my.y.lab = by(data = my.y, INDICES = smpl_table$pop.idx, FUN = mean)


painting.cols = c("#F0E442", "grey60", "#009E73")

max.scaf.pos = round((fai$V2[fai$V1 == my.scaffold]-5e5)/1e6)

png(paste0(res.dir, "/ChromPaint_", my.scaffold, "_",
           length(aq.alleles.aim.noNAmales), "AIMs_", thresh*100, "percDiff_noFinn_pch15_col.png"),
    width = 10, height = 4, units = "in", res = 800)

par(mar = c(4, 6, 1, 1))

plot(real.pos.aim.noNAmales[c(1, length(real.pos.aim.noNAmales))],
     y = c(0, max(my.y)), col = "white", axes = F, xlab = "Position (Mb)", ylab = NA)

# x-axis
axis(1, at = seq(0, max.scaf.pos, by = 1)*1e6, labels = seq(0, max.scaf.pos, by = 1))

# populations
mtext(side = 2, at = my.y.lab[c(3, 4, 6, 7)], #c(20.5, 42, 62.5, 80.5),
      text = c("Bunkkeri", "LangholmenW", "LangholmenR", "Pikkala"), las = 2)
mtext(side = 2, at = my.y.lab[c(2, 8)], #c(20.5, 42, 62.5, 80.5),
      text = c(expression(Finnish~italic(F.~aquilonia)), expression(Finnish~italic(F.~polyctena))), las = 2, cex = .65)
mtext(side = 2, at = my.y.lab[c(1, 9)], text = c(expression(italic(F.~aquilonia)), expression(italic(F.~polyctena))), las = 2)

# data
for (i in 1:ncol(paint.dat)){
  points(real.pos.aim.noNAmales, rep(my.y[i], length(real.pos.aim.noNAmales)),
         col = painting.cols[paint.dat[,i]+1], pch = 15, cex = .4)
}

graphics.off()


#
##
### This is the end. ----
