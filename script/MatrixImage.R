#==============================================================================
# MatrixImage.R : Convert matrix to image
#==============================================================================

# data filename
filenames = dir(path = ".", pattern = "^PI.0.0.mat$")
print(filenames)
for(datafilename in filenames) {
    #______________________________________________________________________________
    # load data from file
    data = read.table(datafilename)

    #______________________________________________________________________________
    # plot matrix as image to default output: screen
    par(pty = 's')
    image(as.matrix(data), axes = T, main = datafilename, col = gray(0:256 / 256))

    #______________________________________________________________________________
	# set png as output target
	png(
        filename = paste(datafilename, ".png", sep = ""),
        width = 1000, height = 1000, units = "px")

    # set postscript as output target
    #postscript(
    #          file = paste(datafilename, ".ps", sep = ""),
    #          width = 8,
    #          height = 8
    #)

	# print output
    par(pty = 's')
    image(as.matrix(data), axes = F, main = datafilename, col = gray(0:256 / 256))
    dev.off()
}
