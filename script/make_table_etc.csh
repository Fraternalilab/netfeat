#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATRUN="$PROJECT/netFeat/network"
set NETFEATTABLE="$PROJECT/netFeat/table"

# clean target
rm $NETFEATTABLE/table_nodes.dat
rm $NETFEATTABLE/table_interactions.dat
rm $NETFEATTABLE/table_k_max.dat 
rm $NETFEATTABLE/table_k_av.dat
rm $NETFEATTABLE/table_k_corr.dat
rm $NETFEATTABLE/table_P_conn.dat

# collect data lines
cd $NETFEATRUN
foreach file (`ls $NETFEATRUN`)
	echo "$file"
	cat $NETFEATRUN/$file/etc.0.dat | grep 'nodes' >> $NETFEATTABLE/table_nodes.dat
	cat $NETFEATRUN/$file/etc.0.dat | grep 'interactions' >> $NETFEATTABLE/table_interactions.dat
	cat $NETFEATRUN/$file/etc.0.dat | grep 'k_max' >> $NETFEATTABLE/table_k_max.dat
	cat $NETFEATRUN/$file/etc.0.dat | grep '<k>' >> $NETFEATTABLE/table_k_av.dat
	cat $NETFEATRUN/$file/etc.0.dat | grep '<kk>' >> $NETFEATTABLE/table_k_corr.dat
	cat $NETFEATRUN/$file/etc.0.dat | grep 'P_conn' >> $NETFEATTABLE/table_P_conn.dat
end
