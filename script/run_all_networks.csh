#! /bin/csh
##! /usr/bin/csh
# run netfeat and R-scripts over all networks
# paths
set PROJECT='~/jk.science/project/ff/netFeat'
set NETWORK="$PROJECT/network"
set NETFEATHOME="$PROJECT/netFeat"
set NETFEATPROG="$PROJECT/netFeat/program/src/netfeat"
set RSCRIPT="$PROJECT/netFeat/script/MatrixImage.R"

# network sub-directories
set files = `ls -1 $NETWORK`

foreach file ($files)
	# create network directory
	#mkdir $NETFEATHOME/network/$file
	# go to directory
	cd $NETFEATHOME/network/$file
	# link network data
	rm -f *
	ln -s $NETWORK/$file/$file.dat
	# run netfeat
	#$NETFEATPROG --matIn $file.dat
	#valgrind  --leak-check=full $NETFEATPROG --matIn $file.dat
	$NETFEATPROG --matIn $file.dat
	# create matrix image
	#R CMD BATCH $RSCRIPT
end
