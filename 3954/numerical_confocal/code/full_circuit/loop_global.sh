# for i in {0..30000..1000}
#   do
#      echo "Welcome $i times"
#  done


for i in $(seq 1 1000 10000)
# for i in $(seq 0 1 0)
# for i in $(seq 1 500 1000)
do
  qsub -v start=$i,end=$(($i+1000)) global.pbs
  # qsub -v start=$i,end=$(($i+1)) global.pbs
  echo "start $i end $(($i+1000))"
done
