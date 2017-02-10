#! /bin/csh
##! /usr/bin/csh
# run netfeat and R-scripts over all networks
# paths
set PROJECT='~/jk.science/project/ff/netFeat'
set NETWORK="$PROJECT/network"
set REFDATA="$PROJECT/reference_data"
set NETFEATPROG="$PROJECT/netFeat/src/netfeat"

# network sub-directories
set files = `ls -1 $NETWORK`

foreach file ($files)
	# go to directory
	cd $NETWORK/$file
	# run netfeat
	#$NETFEATPROG --intsList $REFDATA/$file/intslist.dat --protList $REFDATA/$file/protlist.dat
	diff $file.dat ints.mat.out.0 > diff.log
end
