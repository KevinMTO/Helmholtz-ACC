module purge 

export MODULEPATH=~training/NVIDIA_HPC_SDK_early_access/install/test/modulefiles:$MODULEPATH

module load nvcompilers/2020
module load openmpi/3.1.5/2020

nvidia-smi

array=(1 2 3 4 5 6 7 8)

echo "RUNNING IN GPU MODE"

for i in "${array[@]}"
do
    echo "${i} MPI executed_X"
    echo
	
    name="${i}_multi_mono_X_strong.out"
    output="../WEAK/acc/X/${name}"


    touch $output
    for i in {1..7}
    do
        mpirun -np $i  ../topoGPU_X.exe < ../input >> $output
    done
    
done

export ACC_NUM_CORES=4 
echo "RUNNING IN MULTICORE MODE"
for i in "${array[@]}"
do
    echo "${i} MPI executed_X"
    echo
	
    name="${i}_multi_mono_X_strong.out"
    output="../WEAK/mix/X/${name}"


    touch $output
    for i in {1..7}
    do
        mpirun  -np $i  ../topoMULTICORE_X.exe < ../input >> $output
    done
    
done

echo "RUNNING IN HOST MODE"
for i in "${array[@]}"
do
    echo "${i} MPI executed_X"
    echo
	
    name="${i}_multi_mono_X_strong.out"
    output="../WEAK/pure/X/${name}"


    touch $output
    for i in {1..7}
    do
        mpirun --map-by socket --bind-to core -np $i  ../topoHOST_X.exe < ../input >> $output
    done
    
done