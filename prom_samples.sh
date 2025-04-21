samples=$1
L=$2
Nrun=$3

echo "samples="$samples

for((i=1;i<=$samples;i++))
do 
    echo $i; 
    ./activeinterface $L $Nrun $i; 
    cp inst_sofq.dat "inst_sofq_"$i"_.dat"; 
    cp cm.dat "cm_"$i"_.dat";
done

