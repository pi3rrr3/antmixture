# LAI_loter_vs_genomic_features_per_window.R

# Loter coding: 0 = aquilonia, 1 = polyctena

options(stringsAsFactors = F)
setwd("~/Documents/WORK/speciant/reseq/analysis/local_ancestry/loter")

# Libs
library(GenomicRanges)
library(genomeIntervals)
library(HelloRanges)
library(scales)
assignToRegion = function(bed, pos) {
  #assigns positions (like SNPs) to non-overlapping regions defined in a data frame in bed-style format (scaffold, start, end) and SNPs is a data frame specifying scaffold and position
  pos[,2] = as.numeric(pos[,2])  #makes sure the positions are numbers
  bed[,1] = as.vector(bed[,1])
  pos[,1] = as.vector(pos[,1])    #make scaffold names as vectors or else we cannot compare the names
  #we need to order the bed and pos by scaffold and position
  bed$index = 1:nrow(bed)
  bed =bed[order(bed[,1],bed[,2]),]
  pos$index = 1:nrow(pos)    #to keep track of the original order
  pos = pos[order(pos[,1], pos[,2]),]
  bed[,1] -> bedScaff
  bed[,2] -> start
  bed[,3] -> end
  scaffold = pos[,1]
  posi = pos[,2]
  regIndex = bed$index
  print("assigning to regions...")
  region = rep(NA,length(scaffold))
  j=1
  i=1
  while (i <=length(posi) & j<=nrow(bed)) {
    if (scaffold[i] < bedScaff[j]) {
      i=i+1
      next
    }
    if (scaffold[i] > bedScaff[j] | posi[i]>end[j]) {
      j=j+1
    } else {
      if (posi[i] >= start[j]) {
        region[i]=regIndex[j]
      }
      i=i+1
    }
  }
  region = region[order(pos$index)]
  return(region)
}

# For plotting purposes
my.pops = c("bun", "pik", "lanw", "lanr")
my.full.pops = c("Bunkkeri", "Pikkala", "LångholmenW", "LångholmenR")
my.cols = c("#CC79A7", "#E69F00", "#0072B2", "#D55E00")





###
#### Read LAI results ----
###

lai = read.table("positions_avg_pop_ancestries_nonFin.tab.gz", h = T)





###
#### Gather recombination estimates over 20 kb windows ----
###

# iSMC outputs one bedgraph file per metric, the function below creates a single file with various metrics

iSMCoutput = function(size, my.path, metric){
  tmp = read.table(paste0(my.path, metric[1], ".", size, ".bedgraph"), h = T)
  names(tmp)[4:15] = paste0(metric[1], ".", names(tmp)[4:15])

  for (metricj in metric[-1]){
    tmp2 = read.table(paste0(my.path, metricj, ".", size, ".bedgraph"), h = T)
    names(tmp2)[4:15] = paste0(metricj, ".", names(tmp2)[4:15])
    if(all(tmp$chromStart == tmp2$chromStart)){
      tmp = cbind(tmp, tmp2[,-c(1:3)])
    }else{print("Coordinates not matching...")}
    
  }

  tmp = tmp[, -grep(names(tmp), pattern = "rho"),]
  tmp = tmp[, -which(names(tmp) == "joint_decode")]
  tmp$missing.prop = rowMeans(tmp[,grep(names(tmp), pattern = "missing.prop.")])
  tmp = tmp[, -grep(names(tmp), pattern = "missing.prop.")]
  tmp$diversity.Faq = rowMeans(tmp[,5:10])
  tmp$diversity.Fpol = rowMeans(tmp[,11:16])
  tmp = tmp[, -grep(names(tmp), perl = T, pattern = "diversity.+w")]

  return(tmp)
}

window.res.20k = iSMCoutput(size = "20kb",
                            my.path = "/Users/pnouhaud/Documents/WORK/speciant/reseq/analysis/recomb/final/all_inds.",
                            metric = c("rho", "diversity", "missing.prop"))





###
#### Collect LAI results per window ----
###


# Averages loter results over 20 kb windows (0 = aquilonia, 1 = polyctena)

averageLAI = function(tmp){
  
  lai.key = assignToRegion(tmp[,1:3], lai[,1:2])
  win.key = sort(unique(lai.key))
  
  tmp$n_snp = NA
  tmp$n_snp[win.key] = table(lai.key)
  
  tmp$mean_ancestry_bun = NA
  tmp$mean_ancestry_bun[win.key] = unlist((by(data = lai$ancestry_bun, INDICES = lai.key, FUN = mean)))
  
  tmp$mean_ancestry_pik = NA
  tmp$mean_ancestry_pik[win.key] = unlist((by(data = lai$ancestry_pik, INDICES = lai.key, FUN = mean)))
  
  tmp$mean_ancestry_lanw = NA
  tmp$mean_ancestry_lanw[win.key] = unlist((by(data = lai$ancestry_lanw, INDICES = lai.key, FUN = mean)))
  
  tmp$mean_ancestry_lanr = NA
  tmp$mean_ancestry_lanr[win.key] = unlist((by(data = lai$ancestry_lanr, INDICES = lai.key, FUN = mean)))
  
  
  print("# SNPs per window")
  print(summary(tmp$n_snp))
  
  print(paste0("Fraction of SNPs assigned to windows: ", round(sum(tmp$n_snp, na.rm = T)/nrow(lai)*100, 5), "%"))
  
  hist(tmp$n_snp, breaks = 50)
  rm(win.key, lai.key)
  
  return(tmp)
  
}

window.res.20k = averageLAI(window.res.20k)


# Manual check
bob = window.res.100k[1974,] ; bob
colMeans(lai[which(lai$CHR == bob$chrom & lai$POS >= bob$chromStart & lai$POS <= bob$chromEnd),-c(1:2)])





###
#### Collect CDS fraction per window ----
###

# Write bed with window coordinates, requires bedtools and a gff containing only one transcript per gene

CDSfraction = function(input, wind.size){
  
  # Write bed file
  options(scipen=10)
  write.table(file = "temp.bed", x = cbind(input[,1], input[,2:3]-1),
              col.names = F, row.names = F, quote = F, sep = "\t")
  options(scipen=0)
  
  # Get overlap with bedtools intersect
  system(paste0("bedtools intersect -a temp.bed -b Formica_hybrid_v1_only_longest_CDS.gff3 -wao | cut -f1-3,13 > temp_res.tab"))
  
  # Read results
  temp = read.table("temp_res.tab", h = F)
  names(temp) = c("CHR", "wind.st", "wind.end", "overlap")
  
  # Add overlaps within the same window
  tempkey = paste(temp$CHR, temp$wind.st+1, temp$wind.end+1, sep = "_")
  res = by(data = temp$overlap, INDICES = tempkey, FUN = sum)
  res2 = res[match(unique(tempkey), names(res))]/wind.size # "by" messes with the order of windows, since tempkey is treated as a factor
  summary((paste(input$chrom, input$chromStart, input$chromEnd, sep = "_")) == names(res2))
  input$cds_fraction = res2
  rm(temp, tempkey, res, res2)
  system(paste0("rm temp.bed temp_res.tab"))
  hist(input$cds_fraction)
  return(input)
  
  
}

window.res.20k = CDSfraction(window.res.20k, 20e3)





###
#### Save & load tables ----
###

write.table(window.res.20k, "windows_20kb_rho_ancestry_CDSfract.tsv",
             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

window.res.20k = read.table("/Users/pnouhaud/Documents/WORK/speciant/reseq/analysis/local_ancestry/loter/windows_20kb_rho_ancestry_CDSfract.tsv", header = TRUE)





###
#### Window to IGV ----
###

# if needed, to check alignments etc with IGV

tmp = cbind(window.res.20k[,1], round(rowMeans(window.res.20k[,2:3]))-1, round(rowMeans(window.res.20k[,2:3]))+1, round(window.res.20k[,10], 5))
head(tmp)
cat("track type=bedGraph\n", file="../../../../pacbio/parents/SVs/LocAnc_Pikkala_20kb.bedgraph")
write.table(tmp[complete.cases(tmp),], "../../../../pacbio/parents/SVs/LocAnc_Pikkala_20kb.bedgraph",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append=TRUE)





###
#### Correlation between local ancestry landscapes ----
###

# Check whether average local ancestries are correlated across hybrid populations

png(filename = paste0("plots/correlation_pops_wind20kb.png"),
    w = 12, h = 12, units = "in", res = 200)
par(mfrow = c(4, 4), mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1), cex = 1.5)
for(i in 1:4){
  for(j in 1:4){
    
    if(j >= i){
      
      if (i == j){
        plot(1, 1, col = "white", axes = F, xlab = NA, ylab = NA)
        text(1, 1, my.pops[i], cex = 1.5)
      }else{
        plot(window.res.20k[,8+j], window.res.20k[,8+i],
             xlab = "", ylab = "", pch = 16, col = grey(.1,.1), cex = .5)
        abline(lm(window.res.20k[,8+i] ~ window.res.20k[,8+j]), lwd = 2, col = "orangered")
      }
      
    }else{
        plot(1, 1, col = "white", axes = F, xlab = NA, ylab = NA)
        mycor = cor(window.res.20k[,8+j], window.res.20k[,8+i], method = "spearman", use = "pairwise.complete.obs")
        text(1, 1, round(mycor, 3), cex = 1.2)
      }
  }
}
dev.off()





###
#### Rho quartiles ----
###

nquantile = 5
for (win.size in 20){

  my.dat = get(paste0("window.res.", win.size, "k"))
  
  bob = my.dat[which(my.dat$n_snp >= win.size),]
  bob = within(my.dat, quantile <- as.integer(cut(sample_mean, quantile(sample_mean, probs=0:nquantile/nquantile), include.lowest=TRUE)))
  
  # Folded ancestry
  png(paste0("plots/Rho_quartiles/bp_", win.size, "kb_windows_rhoQuantiles_vs_folded_ancestry.png"), h = 7.5, w = 7.5, units = "in", res = 200)
  par(mar = c(4, 4, 2, 2), oma = c(1, 1, 0, 0), cex = 1.5, mfrow = c(2, 2))
  for(i in 1:4){
    bp = boxplot(abs(bob[,8+i]-.5) ~ bob$quantile, lwd = 1,
                 xlab = "Recombination rate quartile", ylab = paste0("Sorting - ", my.pops[i]))
    for(qt in 1:nquantile){
      points(jitter(rep(qt, sum(bob$quantile == qt)), amount = .35),
             abs(bob[,8+i]-.5)[which(bob$quantile == qt)],
             pch = 16, col = grey(.5, .05), cex = .75)
    }
    points(1:nquantile,
           by(data = abs(bob[,8+i]-.5), INDICES = bob$quantile, FUN = function(x) mean(x, na.rm = TRUE)),
           col = "orangered2", pch = 16, cex = 1.2)
  }
  dev.off()
  
  # Full ancestry
  png(paste0("plots/Rho_quartiles/bp_", win.size, "kb_windows_rhoQuantiles_vs_ancestry.png"), h = 7.5, w = 7.5, units = "in", res = 200)
  par(mar = c(4, 4, 2, 2), oma = c(1, 1, 0, 0), cex = 1.5, mfrow = c(2, 2))
  for(i in 1:4){
    bp = boxplot(bob[,8+i] ~ bob$quantile, lwd = 1,
                 xlab = "Recombination rate quartile", ylab = paste0("Ancestry - ", my.pops[i]))
    for(qt in 1:nquantile){
      points(jitter(rep(qt, sum(bob$quantile == qt)), amount = .35),
             bob[,8+i][which(bob$quantile == qt)],
             pch = 16, col = grey(.5, .05), cex = .75)
    }
    points(1:nquantile,
           by(data = bob[,8+i], INDICES = bob$quantile, FUN = function(x) mean(x, na.rm = TRUE)),
           col = "orangered2", pch = 16, cex = 1.2)
  }
  dev.off()
  
}






###
#### CDS quartiles ----
###

nquantile = 4
for (win.size in c(20, 50, 100)){
  
  my.dat = get(paste0("window.res.", win.size, "k"))
  
  bob = my.dat[which(my.dat$n_snp >= win.size),]
  bob = within(bob, quantile <- as.integer(cut(cds_fraction, quantile(cds_fraction, probs=0:nquantile/nquantile), include.lowest=TRUE)))
  quantile(bob$cds_fraction, 0:nquantile/nquantile)
  
  # Folded ancestry
  png(paste0("plots/CDS_quartiles/bp_", win.size, "kb_windows_CDSQuantiles_vs_folded_ancestry.png"), h = 7.5, w = 7.5, units = "in", res = 200)
  par(mar = c(4, 4, 2, 2), oma = c(1, 1, 0, 0), cex = 1.5, mfrow = c(2, 2))
  for(i in 1:4){
    bp = boxplot(abs(bob[,8+i]-.5) ~ bob$quantile, lwd = 1,
                 xlab = "CDS fraction quartile", ylab = paste0("Sorting - ", my.pops[i]))
    for(qt in 1:nquantile){
      points(jitter(rep(qt, sum(bob$quantile == qt)), amount = .35),
             abs(bob[,8+i]-.5)[which(bob$quantile == qt)],
             pch = 16, col = grey(.5, .05), cex = .75)
    }
    points(1:nquantile,
           by(data = abs(bob[,8+i]-.5), INDICES = bob$quantile, FUN = function(x) mean(x, na.rm = TRUE)),
           col = "orangered2", pch = 16, cex = 1.2)
  }
  dev.off()
  
  # Full ancestry
  png(paste0("plots/CDS_quartiles/bp_", win.size, "kb_windows_CDSQuantiles_vs_ancestry.png"), h = 7.5, w = 7.5, units = "in", res = 200)
  par(mar = c(4, 4, 2, 2), oma = c(1, 1, 0, 0), cex = 1.5, mfrow = c(2, 2))
  for(i in 1:4){
    bp = boxplot(bob[,8+i] ~ bob$quantile, lwd = 1,
                 xlab = "CDS fraction quartile", ylab = paste0("Ancestry - ", my.pops[i]))
    for(qt in 1:nquantile){
      points(jitter(rep(qt, sum(bob$quantile == qt)), amount = .35),
             bob[,8+i][which(bob$quantile == qt)],
             pch = 16, col = grey(.5, .05), cex = .75)
    }
    points(1:nquantile,
           by(data = bob[,8+i], INDICES = bob$quantile, FUN = function(x) mean(x, na.rm = TRUE)),
           col = "orangered2", pch = 16, cex = 1.2)
  }
  dev.off()
  
}




#
##
### This is the end.
