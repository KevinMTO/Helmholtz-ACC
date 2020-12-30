    #!bin/bash
    export ACC_NUM_CORES=2


    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 1 ../topoGPU_alongY.exe < input1
    done >> 1_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 2 ../topoGPU_alongY.exe < input2
    done >> 2_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 3 ../topoGPU_alongY.exe < input3
    done >> 3_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 4 ../topoGPU_alongY.exe < input4
    done >> 4_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 5 ../topoGPU_alongY.exe < input5
    done >> 5_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 6 ../topoGPU_alongY.exe < input6
    done >> 6_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 7 ../topoGPU_alongY.exe < input7
    done >> 7_GPU_alongY.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 8 ../topoGPU_alongY.exe < input8
    done >> 8_GPU_alongY.txt




    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 1 ../topoGPU_alongX.exe < input1
    done >> 1_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 2 ../topoGPU_alongX.exe < input2
    done >> 2_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 3 ../topoGPU_alongX.exe < input3
    done >> 3_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 4 ../topoGPU_alongX.exe < input4
    done >> 4_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 5 ../topoGPU_alongX.exe < input5
    done >> 5_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 6 ../topoGPU_alongX.exe < input6
    done >> 6_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 7 ../topoGPU_alongX.exe < input7
    done >> 7_GPU_alongX.txt

    for j in {1..6}
    do
        mpirun --map-by socket --bind-to core -np 8 ../topoGPU_alongX.exe < input8
    done >> 8_GPU_alongX.txt