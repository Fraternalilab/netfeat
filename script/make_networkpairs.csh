#! /usr/bin/csh
##! /bin/csh

# create paired network names/directories
# paths
#set HOME=/Users/jkleinj
set HOME=/home/jkleinj
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETWORK="$PROJECT/network"
set NETWORKPAIR="$PROJECT/netFeat/networkpair"

# network sub-directories
set files = `ls -1 $NETWORK`

# create paired network names/directories
set i = 0
foreach file0 ($files)
	@ i ++
	set j = 0
	foreach file1 ($files)
		@ j ++
		if ($j < $i) then
			echo "mkdir $NETWORKPAIR/$file0-$file1"
			mkdir "$NETWORKPAIR/$file0-$file1"
		endif
	end
end

