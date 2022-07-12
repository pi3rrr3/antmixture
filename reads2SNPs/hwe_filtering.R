# hwe_filtering.R

# Get the filename to read the hwe test
filename = "FEMALES_hwe" # file name of reference table

# Read the HW results
hwvalues = read.table(paste(filename,".hwe",sep=""), header=T, stringsAsFactors = F)

# Get p-value distributions
print(paste("he excess", sum(hwvalues$P_HET_EXCESS<0.01)))
print(paste("he deficit", sum(hwvalues$P_HET_DEFICIT<0.01)))
print(paste("he overall", sum(hwvalues$P_HWE<0.01)))

# Remove only the sites with excess of heterozygotes
idx.remove = which(hwvalues$P_HET_EXCESS<0.01)

# Create BED file with the sites that fail the HW test for excess of heterozygotes
position = hwvalues[idx.remove, 2]
bedmatrix = matrix(c(hwvalues[indextoremove, 1], format(position-1, scientific = F), format(position, scientific = F)),
	    	   ncol = 3, byrow = F)
write("#Sites with heterozygosity excess", file = paste("nohwe_excess_", filename, ".bed", sep=""))
write(t(bedmatrix), file = paste("nohwe_excess_", filename, ".bed", sep=""), ncolumns = 3, append = T)