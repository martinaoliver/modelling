L=$1
dx=$2
T=$3
dt=$4
python colonyMaskCreation.py $L $dx $T $dt
python run_numerical_single.py $L $dx $T $dt
