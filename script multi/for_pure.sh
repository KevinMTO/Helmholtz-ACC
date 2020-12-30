cd#!/bin/sh

export ACC_NUM_CORES=8

for j in {1..6}
do
    mpirun -np 1 --map-by ppr:1:node:pe=8 --bind-to core --report-bindings ../topohost_2D.exe < ../compiled/input1
done >> ./1_MPI8cores2Dpure.txt

export ACC_NUM_CORES=4

for j in {1..6}
do
    mpirun -np 2 --map-by ppr:1:socket:pe=4 --bind-to core --report-bindings ../topohost_2D.exe < ../compiled/input2
done >> ./2_MPI4cores2Dpure.txt

export ACC_NUM_CORES=2

for j in {1..6}
do
    mpirun -np 4 --map-by ppr:2:socket:pe=2 --bind-to core --report-bindings ../topohost_2D.exe < ../compiled/input4
done >> ./4_MPI2cores2Dpure.txt

export ACC_NUM_CORES=1

for j in {1..6}
do
    mpirun -np 8 --map-by ppr:4:socket:pe=1 --bind-to core --report-bindings ../topohost_2D.exe < ../compiled/input8
done >> ./8_MPI1cores2Dpure.txt









export ACC_NUM_CORES=8

for j in {1..6}
do
    mpirun -np 1 --map-by ppr:1:node:pe=8 --bind-to core --report-bindings ../topoHOST_alongX.exe < ../compiled/input1
done >> ./1_MPI8coresXpure.txt

export ACC_NUM_CORES=4

for j in {1..6}
do
    mpirun -np 2 --map-by ppr:1:socket:pe=4 --bind-to core --report-bindings ../topoHOST_alongX.exe < ../compiled/input2
done >> ./2_MPI4coresXpure.txt

export ACC_NUM_CORES=2

for j in {1..6}
do
    mpirun -np 4 --map-by ppr:2:socket:pe=2 --bind-to core --report-bindings ../topoHOST_alongX.exe < ../compiled/input4
done >> ./4_MPI2coresXpure.txt

export ACC_NUM_CORES=1

for j in {1..6}
do
    mpirun -np 8 --map-by ppr:4:socket:pe=1 --bind-to core --report-bindings ../topoHOST_alongX.exe < ../compiled/input8
done >> ./8_MPI1coresXpure.txt









export ACC_NUM_CORES=8

for j in {1..6}
do
    mpirun -np 1 --map-by ppr:1:node:pe=8 --bind-to core --report-bindings ../topoHOST_alongY.exe < ../compiled/input1
done >> ./1_MPI8coresYpure.txt

export ACC_NUM_CORES=4

for j in {1..6}
do
    mpirun -np 2 --map-by ppr:1:socket:pe=4 --bind-to core --report-bindings ../topoHOST_alongY.exe < ../compiled/input2
done >> ./2_MPI4coresYpure.txt

export ACC_NUM_CORES=2

for j in {1..6}
do
    mpirun -np 4 --map-by ppr:2:socket:pe=2 --bind-to core --report-bindings ../topoHOST_alongY.exe < ../compiled/input4
done >> ./4_MPI2coresYpure.txt

export ACC_NUM_CORES=1

for j in {1..6}
do
    mpirun -np 8 --map-by ppr:4:socket:pe=1 --bind-to core --report-bindings ../topoHOST_alongY.exe < ../compiled/input8
done >> ./8_MPI1coresYpure.txt