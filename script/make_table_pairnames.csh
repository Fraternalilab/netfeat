#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETWORKPAIR="$PROJECT/netFeat/networkpair"

# network sub-directories
set files = `ls -1 $NETWORKPAIR`

set i = 0
foreach file ($files)
	echo "$file"
end

