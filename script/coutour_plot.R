# data filename and starting res
filenames = dir(path = ".", pattern = "^PI\.*mat$")
print(filenames)


for(datafilename in filenames) {
 data = as.matrix(read.table(datafilename)) 	 	    
 num_rows=nrow(data)
 kmax=40
 if(num_rows <40){
  kmax=num_rows
 }
 a=1
 b=seq(5,kmax,by=5)
 a=append(a,b)
 a=append(a,kmax)
 # print postscript file
 pdf(
  file = paste(datafilename, ".pdf", sep = ""),width = 7,height = 6.4
 )
 filled.contour(x=seq(1,kmax),y=seq(1,kmax),data[1:kmax,1:kmax],levels=seq(0,1.5,by=0.02),color.palette=colorRampPalette(c("black","blue","light blue","yellow","red","white")),key.axes= axis(4, seq(0, 1.5, by = 0.1)),key.title = title(main="Pi"),xlab="k",ylab="k`",asp=1,frame.plot=FALSE,plot.axes = { axis(1, at=a);axis(2, at=a)})
 dev.off()
}
