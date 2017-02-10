#==============================================================================
# MatrixImage.R : Convert matrix to image
# (C) 2008 Jens Kleinjung
#==============================================================================

# data filename and starting res
filenames = dir(path = ".", pattern = "^kernel\.*dat$")
print(filenames)
for(datafilename in filenames) {
    #______________________________________________________________________________
    # map matrix as image 
    # load data from file
    data = read.table(datafilename)

    #______________________________________________________________________________
    # plot map
    par(pty = 's')
    #image(as.matrix(data), axes = F, main = datafilename)
    image(as.matrix(data), axes = T, main = datafilename, col = gray(0:256 / 256))

    #______________________________________________________________________________
    # print postscript file
    postscript(
              file = paste(datafilename, ".ps", sep = ""),
              width = 8,
              height = 8
    )
    par(pty = 's')
    image(as.matrix(data), axes = F, main = datafilename)
    dev.off()
}
