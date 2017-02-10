#! /usr/bin/csh
##! /bin/csh

# create paired network names/directories
# paths
#set HOME=/Users/jkleinj
set HOME=/home/jkleinj
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETWORK="$PROJECT/network"
set NETWORKDATA="../../network"
set NETWORKPAIR="$PROJECT/netFeat/networkpair"

# network sub-directories
set files = `ls -1 $NETWORK`

# link result data to paired network names/directories
set i = 0
foreach file0 ($files)
	@ i ++
	set j = 0
	foreach file1 ($files)
		@ j ++
		if ($j < $i) then
			echo "$NETWORKPAIR/$file0-$file1"
			cd "$NETWORKPAIR/$file0-$file1"
			ln -fs $NETWORKDATA/$file0/$file0.dat
			ln -fs $NETWORKDATA/$file1/$file1.dat
		endif
	end
end

