# loter_preprocessing.R

options(stringsAsFactors = F)
library(data.table)
setwd("~/Documents/WORK/speciant/reseq/analysis/local_ancestry/loter")



# For plotting purposes
my.pops = c("bun", "pik", "lanw", "lanr")
my.full.pops = c("Bunkkeri", "Pikkala", "LångholmenW", "LångholmenR")
my.cols = c("#CC79A7", "#E69F00", "#0072B2", "#D55E00", "#009E73", "#F0E442")
names(my.cols) = c("bun", "pik", "lanw", "lanr", "Faq", "Fpol")




###
#### Read data ----
###

# 0 = aquilonia, 1 = polyctena

positions = fread("positions.tab.gz", h = F, col.names = c("CHR", "POS"))
bun = fread("loter_nonFin.bun.females.txt.gz", h = F)
pik = fread("loter_nonFin.pik.females.txt.gz", h = F)
lanw = fread("loter_nonFin.lanw.females.txt.gz", h = F)
lanr = fread("loter_nonFin.lanr.females.txt.gz", h = F)

res = read.table("tract_lengths_loter_nonFin.tsv.gz", h=T)

par(mfrow = c(4, 1))
hist(res$length[res$population == "bun"], breaks = 1000, xlim = c(0, 5e5), col = "black")
hist(res$length[res$population == "pik"], breaks = 1000, xlim = c(0, 5e5), col = "black")
hist(res$length[res$population == "lanw"], breaks = 1000, xlim = c(0, 5e5), col = "black")
hist(res$length[res$population == "lanr"], breaks = 1000, xlim = c(0, 5e5), col = "black")




###
#### Compute tract lengths ----
###

scaf = unique(positions$CHR)[1]
res = NULL
for (popname in c("bun", "pik", "lanw", "lanr")){
  
  print(paste0("###### Population ", popname))
  pop = get(popname)
  
  for (scaf in unique(positions$CHR)){
    
    print(paste0("### ", scaf))
    pop.scaf = pop[, which(positions$CHR == scaf), with = F]
    
    for(i in 1:nrow(pop.scaf)){
      
      print(paste0("Haplotype ", i))
      
      start = which(sequence(rle(pop.scaf[i])$lengths) == 1)
      end = start + rle(pop.scaf[i])$lengths - 1
      if(pop.scaf[i, 1, with=F] == 1){
        ancestry = rep(c(1, 0), length(start))[1:length(start)]
      }else{
        ancestry = rep(c(0, 1), length(start))[1:length(start)]
      }
      
      tmp = rbind(data.frame(population = popname,
                             scaffold = scaf,
                             haplotype = i,
                             start = positions$POS[which(positions$CHR == scaf)][start],
                             end = positions$POS[which(positions$CHR == scaf)][end],
                             ancestry,
                             length = positions$POS[which(positions$CHR == scaf)][end] - positions$POS[which(positions$CHR == scaf)][start]))
      
      res = rbind(res, tmp)
      
      hist(tmp$length, main = paste(popname, scaf, i))
      
    }
    
    rm(pop.scaf, tmp, start, end) ; gc()
    
  }
  
  rm(pop) ; gc()
  
}

write.table(res, "tract_lengths_loter_nonFin.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


# Tract length sizes

png(filename = paste0("plots/hist_tract_length_per_haplo_min1kb.png"),
    w = 12, h = 12, units = "in", res = 200)
par(mfrow = c(1, 1), cex = 1.5)
plot(1, 1, ylim = c(0, 120), xlim = c(0, 3e5), xlab = "Tract length", ylab = "Count per haploid genome", col = "white")
for (pop in unique(res$population)){
  nhaplo = length(unique(res$haplotype[res$population == pop]))
  tmp = hist(res$length[res$population == pop & res$length > 1e3], breaks = 1000, plot = F)
  lines(tmp$mids, tmp$counts/nhaplo, col = my.cols[which(pop == unique(res$population))], lwd = 4)
  rm(tmp)
}
legend("topright", legend = my.full.pops, lwd = 4, col = my.cols)
graphics.off()


png(filename = paste0("plots/hist_polarized_tract_length_per_haplo_min1kb.png"),
    w = 12, h = 12, units = "in", res = 200)
par(mfrow = c(2, 2), cex = 1.5)
for (pop in unique(res$population)){
  nhaplo = length(unique(res$haplotype[res$population == pop]))
  plot(1, 1, ylim = c(0, 150), xlim = c(0, 5e5), xlab = "Tract length", ylab = "Counts per haploid genome", col = "white", main = my.full.pops[which(pop == unique(res$population))])
  tmp = hist(res$length[res$population == pop
                                      & res$length > 1e3
                                      & res$ancestry == 0], breaks = 500, plot = F)
  lines(tmp$mids, tmp$counts/nhaplo, col = my.cols2[5], lwd = 4, lty = 2)
  tmp2 = hist(res$length[res$population == pop
                                      & res$length > 1e3
                                      & res$ancestry == 1], breaks = 500, plot = F)
  lines(tmp2$mids, tmp2$counts/nhaplo, col = my.cols2[6], lwd = 4, lty = 2)
  rm(tmp, tmp2)
}
dev.off()




###
#### Compute diploid genotypes ----
###

haplo2geno = data.frame(geno = 1:10, haplo1 = seq(1, 19, by = 2), haplo2 = seq(1, 19, by = 2)+1)

geno = matrix(nrow = nrow(positions), ncol = 39, data = NA)

j = 1 # counts inds in geno

for (popname in c("bun", "pik", "lanw", "lanr")){
  
  print(paste0("###### Population ", popname))
  pop = get(popname)
  nind = nrow(pop)/2
  for(i in 1:nind){
    geno[,j] = colSums(pop[c(haplo2geno$haplo1[i], haplo2geno$haplo2[i]),])
    j = j+1
  }

}

write.table(geno, "genotypes_from_loter.tsv", sep = "\t", col.names = F, row.names = F, quote = F)




###
#### Compute average LA per haplotype ----
###

LAperhaplo = matrix(ncol = 20, nrow = 4, data = NA)
for (popname in c("bun", "pik", "lanw", "lanr")){
  pop = get(popname)
  LAperhaplo[which(popname == c("bun", "pik", "lanw", "lanr")), 1:nrow(pop)] = rowSums(pop == 1)/ncol(pop)
  }
rownames(LAperhaplo) = c("bun", "pik", "lanw", "lanr")

# Plot

png("plots/average_LA_per_haplotype_nonFin.png", h = 7.5, w = 7.5, units = "in", res = 200)
par(mar = c(4, 4, 1, 1), cex = 1.5, mfrow = c(4, 1))

my.x = sort(c(seq(1, 10, by = 1)-.1, seq(1, 10, by = 1)+.1))
for(i in 1:4){
  nhaplo = sum(!is.na(LAperhaplo[i,]))
  plot(my.x[1:nhaplo], LAperhaplo[i, 1:nhaplo],
       col = "white", type = "p", pch = 16,
       ylim = c(.2, .8), xlim = c(.5, 10),
       axes = F, xlab = "", ylab = "Local ancestry")
  axis(2, at = seq(.2, .8, by = .1), las = 2) ; axis(1, at = c(1:(nhaplo/2)))
  abline(h = seq(0, 1, by = .1), lwd = .2, col = "grey")
  abline(h = .5, lwd = .5, col = "grey")
  points(my.x[1:nhaplo], LAperhaplo[i, 1:nhaplo], cex = 1.5, pch = 16)
  text(1, .7, c("bun", "pik", "lanw", "lanr")[i], cex = 1.5)
}

dev.off()


###
#### Compute average LA per SNP and population ----
###

positions = as.data.frame(positions)
positions$ancestry_bun = NA
positions$ancestry_pik = NA
positions$ancestry_lanw = NA
positions$ancestry_lanr = NA

for (popname in c("bun", "pik", "lanw", "lanr")){
  pop = get(popname)
  tmp = colMeans(pop)
  positions[, grep(pattern = popname, names(positions))] = tmp
  rm(pop)
}

# This file will be used for the subsequent window-based steps
write.table(positions, "positions_avg_pop_ancestries_nonFin.tab", sep = "\t", col.names = T, row.names = F, quote = F)




###
#### Plotting ----
###


my.pops = c("bun", "pik", "lanw", "lanr")
my.full.pops = c("Bunkkeri", "Pikkala", "Langholmen-W", "Langholmen-R")
my.cols = c("#CC79A7", "#E69F00", "#0072B2", "#D55E00")
my.parental.cols = c("#009E73", "#F0E442") # Faq, Fpol

# Chromosome painting with loter

for(scaf in unique(positions$CHR)){
  
  idx.scaf = which(positions$CHR == scaf)
  my.pos = positions$POS[idx.scaf]
  
  print(scaf)
  
  png(file = paste0("plots/nonFin/loter_females_", scaf, "_nonFin.png"), width = 12, height = 6, units = "in", res = 150)
  
  par(mar = c(4, 5, 2, 2))
  
  plot(c(0, max(positions$POS[idx.scaf])), c(0, 93), col = "white", axes = F, xlab = paste0(scaf, " (Mb)"), ylab = NA)
  axis(1, at = seq(0, max(my.pos), by = 1e6), labels = seq(0, max(my.pos), by = 1e6)/1e6)
  axis(2, at = c(10, 35, 60, 84), labels = c("Bunkkeri", "Pikkala", "LangholmenW", "LangholmenR"), cex.axis = .8, tick = F)
  
  j = 0
  
  for (pop in c("bun", "pik", "lanw", "lanr")){
    print(pop)
    tmp = get(pop)[,idx.scaf, with = F]
    j = j + 1
    for (i in 1:nrow(tmp)){
      col.idx = unlist(tmp[i]+1)
      points(x = my.pos, y = rep(j, length(idx.scaf)), col = c("#009E73", "#F0E442")[col.idx], pch = 15, cex = .6)
      j = j + 1
    } ; j = j + 4
  }
  
  dev.off()
}



#
##
### This is the end.
