data = as.matrix(read.table("PI_corr.0.0.mat"))

z = data
x = (1:nrow(z))
y = (1:ncol(z))

# generate colors for 0->1 in blue and 1->x in red
#cellcol = matrix(rep("#000000",nrow(z)))
cellcol = matrix(rep("#000000",20))
cellcol[z<1] = color.scale(z[z<1], c(0,0.8), 0, c(1,0.2))
cellcol[z>=1] = color.scale(z[z>=1], c(1,0.8), 0, c(0.2,0))

color2D.matplot(z,cellcolors=cellcol,show.legend=TRUE,xlab="network 1",ylab="network 2",main="PI correlation")

