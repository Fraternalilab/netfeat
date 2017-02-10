#! /usr/bin/csh

# paths
set PROJECT='/home/jkleinj/jk.science/project/ff/netFeat'
set NETWORK="$PROJECT/network"

# network sub-directories
set files = `ls -1 $NETWORK`

foreach file ($files)
	cd $NETWORK/$file
	echo "$file"
	cat $file.dat | wc -l
end
