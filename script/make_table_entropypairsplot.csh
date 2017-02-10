#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATRUN="$PROJECT/netFeat/network"
set NETFEATTABLE="$PROJECT/netFeat/table"

# clean target
rm $NETFEATTABLE/table_entropypairs.dat

# collect data lines
cd $NETFEATRUN
foreach file (`ls $NETFEATRUN`)
	echo "$file"
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'k_av'         | cut -d ' ' -f 1,2  >> $NETFEATTABLE/table_entropypairs.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'N_entropy'    | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropypairs.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'p_entropy'    | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropypairs.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'PI_entropy'   | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropypairs.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'S_entropy'    | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropypairs.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'C_complexity' | cut -d ' ' -f 2  >> $NETFEATTABLE/table_entropypairs.dat
end
