#! /usr/bin/csh
##! /bin/csh

# paths
#set HOME="/Users/jkleinj"
set HOME="/home/jkleinj"
set PROJECT="$HOME/jk.science/project/ff/netFeat"
set NETFEATRUN="$PROJECT/netFeat/networkpair"
set NETFEATTABLE="$PROJECT/netFeat/table"

# clean target
rm $NETFEATTABLE/table_D_kmax.dat
rm $NETFEATTABLE/table_D_A_PI.dat
rm $NETFEATTABLE/table_D_A_p.dat
rm $NETFEATTABLE/table_D_B_PI.dat
rm $NETFEATTABLE/table_D_B_p.dat
rm $NETFEATTABLE/table_D.dat

# collect data lines
cd $NETFEATRUN
foreach file (`ls $NETFEATRUN`)
	echo "$file"
	cat $NETFEATRUN/$file/KL_distance.0.dat | grep 'k_max' >> $NETFEATTABLE/table_D_kmax.dat
	cat $NETFEATRUN/$file/KL_distance.0.dat | grep 'D_A_p' >> $NETFEATTABLE/table_D_A_p.dat
	cat $NETFEATRUN/$file/KL_distance.0.dat | grep 'D_A_PI' >> $NETFEATTABLE/table_D_A_PI.dat
	cat $NETFEATRUN/$file/KL_distance.0.dat | grep 'D_B_p' >> $NETFEATTABLE/table_D_B_p.dat
	cat $NETFEATRUN/$file/KL_distance.0.dat | grep 'D_B_PI' >> $NETFEATTABLE/table_D_B_PI.dat
	cat $NETFEATRUN/$file/KL_distance.0.dat | grep 'D ' >> $NETFEATTABLE/table_D.dat
end
