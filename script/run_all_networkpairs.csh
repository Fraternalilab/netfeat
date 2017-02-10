#! /bin/csh
##! /usr/bin/csh
# run netfeat over all network pairs

# pre-processing (run only once)
#make_networkpairs.csh
#link_networkpairs.csh

# paths
set files = `ls -1 $NETWORKPAIR`
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATPROG="$PROJECT/netFeat/src/netfeat"
set NETWORK="$PROJECT/network"
set NETWORKDATA="../../network"
set NETWORKPAIR="$PROJECT/netFeat/networkpair"

# network sub-directories
set files = `ls -1 $NETWORK`

set i = 0
foreach file0 ($files)
	@ i ++
	set j = 0
	foreach file1 ($files)
		@ j ++
		if ($j < $i) then
			echo "$NETWORKPAIR/$file0-$file1"
			cd "$NETWORKPAIR/$file0-$file1"
			#valgrind --leak-check=full $NETFEATPROG --matIn $file0.dat --matIn1 $file1.dat --compare
			$NETFEATPROG --matIn $file0.dat --matIn1 $file1.dat --compare 
		endif
	end
end


