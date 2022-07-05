# for i in {0..30000..1000}
#   do
#      echo "Welcome $i times"
#  done

# var_array = (0.1 0.8 0.06 0.04 0.03 0.01 0.001)

# for i in HKS kSI IKS_real IKS_im
for i in 50 100 200 
# for i in 10 20 30 50 60 
# for i in $(seq 0 1 0)
# for i in var_array 
do
  # qsub -v metric=psprenorm,var=$i run_qsub.pbs
  # qsub -v metric=psprenorm run_qsub.pbs
    # qsub -v tg=$i,cf=1 run_qsub.pbs 
    # qsub -v tg=$i,cf=1.1 run_qsub.pbs 
    # qsub -v tg=$i,cf=2 run_qsub.pbs
    qsub -v tg=$i,cf=1.5 run_qsub.pbs 

done
