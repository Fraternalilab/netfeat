#! /usr/bin/csh
# run netfeat and R-scripts over all networks
# paths
set PROJECT='/home/jkleinj/jk.science/project/ff/netFeat'
set NETWORK="$PROJECT/network"
set NETFEATHOME="$PROJECT/netFeat"
set NETFEATPROG="$PROJECT/netFeat/netfeat"
set RSCRIPT="$PROJECT/netFeat/script/MatrixImage.R"

# network sub-directories
set files = `ls -1 $NETWORK`

foreach file ($files)
	# go to directory
	if (-e $NETFEATHOME/network/$file) then
		rm $NETFEATHOME/network/$file/*
	endif
end
