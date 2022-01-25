# for i in {0..30000..1000}
#   do
#      echo "Welcome $i times"
#  done


for i in $(seq 1 1000 10000)
# for i in $(seq 1 1000 30000)
# for i in $(seq 0 1 0)
do
  # qsub -v start=$i,end=$(($i+1)) global.pbs
  qsub -v start=$i,end=$(($i+1000)) global.pbs
   echo "Welcome $i times and $(($i+1))"
done
