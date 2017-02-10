#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATRUN="$PROJECT/netFeat/network"
set NETFEATTABLE="$PROJECT/netFeat/table"

# clean target
#rm $NETFEATTABLE/table_names.dat

# collect data lines
cd $NETFEATRUN
foreach file (`ls $NETFEATRUN`)
	echo "$file"
end
