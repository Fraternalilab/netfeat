#==============================================================================
# Plot all pairs of entropic values in a matrix of scatter plots
#==============================================================================

dat = read.table("table_entropypairs.formatted.dat", header = T)
pairs(dat)

# postscript plot
outFileName = paste("entropypairs.eps") 
postscript(outFileName, width = 7.0, height = 7.0, horizontal = FALSE, onefile = TRUE, paper = "special")
pairs(dat)
dev.off()

