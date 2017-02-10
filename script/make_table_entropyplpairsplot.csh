#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATRUN="$PROJECT/netFeat/network"
set NETFEATTABLE="$PROJECT/netFeat/table"

# clean target
rm $NETFEATTABLE/table_entropyplpairs.dat

# collect data lines
cd $NETFEATRUN
foreach file (`ls $NETFEATRUN`)
	echo "$file"
	cat $NETFEATRUN/$file/entropytab.0.dat | grep 'k_av'            | cut -d ' ' -f 1,2  >> $NETFEATTABLE/table_entropyplpairs.dat
	cat $NETFEATRUN/$file/entropytab.0.dat | grep 'N_entropy_pl'    | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropyplpairs.dat
	cat $NETFEATRUN/$file/entropytab.0.dat | grep 'p_entropy_pl'    | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropyplpairs.dat
	cat $NETFEATRUN/$file/entropytab.0.dat | grep 'PI_entropy_pl'   | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropyplpairs.dat
	cat $NETFEATRUN/$file/entropytab.0.dat | grep 'S_entropy_pl'    | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropyplpairs.dat
	cat $NETFEATRUN/$file/entropytab.0.dat | grep 'C_complexity_pl' | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropyplpairs.dat
end
