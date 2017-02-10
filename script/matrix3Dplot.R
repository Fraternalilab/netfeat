library(plotrix)
data = as.matrix(read.table("PI_corr.0.1.mat"))
#data = as.matrix(read.table("test.dat"))

z = data
x = (1:nrow(z))
y = (1:ncol(z))

# generate colors for 0->1 in blue and 1->x in red
cellcol = matrix(rep("#000000",64))
cellcol[z<1] = color.scale(z[z<1], c(0,1), 0, c(1,0))
cellcol[z>=1] = color.scale(z[z>=1], c(1,0), 0, c(0,1))

persp(x, y, z, theta = 110, phi = 40, col = cellcol, scale = FALSE,
      ltheta = -120, shade = 0.4, border = NA, box = FALSE)

