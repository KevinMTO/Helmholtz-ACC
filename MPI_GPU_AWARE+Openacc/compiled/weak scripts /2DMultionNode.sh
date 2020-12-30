module purge 

export MODULEPATH=~training/NVIDIA_HPC_SDK_early_access/install/test/modulefiles:$MODULEPATH

module load nvcompilers/2020
module load openmpi/3.1.5/2020

nvidia-smi

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

echo "RUNNING IN GPU MODE"

for i in "${array[@]}"
do
    echo "${i} MPI executed_2D"
    echo
	
    name="${i}_multi_mono_2D_strong.out"
    output="../WEAK/acc/2D/${name}"


    touch $output
    for j in {1..7}
    do
        mpirun -np $i  ../topoGPU_2D.exe < ../"input${i}" >> $output
    done
    
done



echo "RUNNING IN HOST MODE"
for i in "${array[@]}"
do
    echo "${i} MPI executed_2D"
    echo
	
    name="${i}_multi_mono_2D_strong.out"
    output="../WEAK/pure/2D/${name}"


    touch $output
    for j in {1..7}
    do
        mpirun --map-by socket --bind-to core -np $i  ../topoHOST_2D.exe < ../"input${i}" >> $output
    done
    
done


echo "RUNNING IN MULTICORE MODE"
for i in "${array[@]}"
do
    export ACC_NUM_CORES=`expr 20 \/ $i` 
    echo "${i} MPI executed_2D"
    echo "number of threads ${ACC_NUM_CORES}"
    echo
	
    name="${i}_multi_mono_2D_strong.out"
    output="../WEAK/mix/2D/${name}"


    touch $output
    for j in {1..7}
    do
        mpirun  -np $i  ../topoMULTICORE_2D.exe < ../"input${i}" >> $output
    done
    
done