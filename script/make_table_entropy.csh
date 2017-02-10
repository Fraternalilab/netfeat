#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATRUN="$PROJECT/netFeat/network"
set NETFEATTABLE="$PROJECT/netFeat/table"

# clean target
rm $NETFEATTABLE/table_S_0.dat
rm $NETFEATTABLE/table_C_p.dat
rm $NETFEATTABLE/table_C_PI.dat 
rm $NETFEATTABLE/table_S.dat
rm $NETFEATTABLE/table_file_network.dat

# collect data lines
cd $NETFEATRUN
foreach file (`ls $NETFEATRUN`)
	echo "$file"
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'S_0' >> $NETFEATTABLE/table_S_0.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'C_p' >> $NETFEATTABLE/table_C_p.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'C_PI' >> $NETFEATTABLE/table_C_PI.dat
	cat $NETFEATRUN/$file/entropy.0.dat | grep 'S ' >> $NETFEATTABLE/table_S.dat
	echo "$file" >> $NETFEATTABLE/table_file_network.dat 
end
