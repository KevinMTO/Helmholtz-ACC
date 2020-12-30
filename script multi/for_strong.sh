#!/bin/sh

export ACC_NUM_CORES=8

for j in {1..6}
do
    mpirun -np 1 --map-by ppr:1:node:pe=8 --bind-to core --report-bindings ../topomulticore_2D.exe < ../compiled/input8
done >> ./strong/1_MPI8cores2D.txt

export ACC_NUM_CORES=4

for j in {1..6}
do
    mpirun -np 2 --map-by ppr:1:socket:pe=4 --bind-to core --report-bindings ../topomulticore_2D.exe < ../compiled/input8
done >> ./strong/2_MPI4cores2D.txt

export ACC_NUM_CORES=2

for j in {1..6}
do
    mpirun -np 4 --map-by ppr:2:socket:pe=2 --bind-to core --report-bindings ../topomulticore_2D.exe < ../compiled/input8
done >> ./strong/4_MPI2cores2D.txt

export ACC_NUM_CORES=1

for j in {1..6}
do
    mpirun -np 8 --map-by ppr:4:socket:pe=1 --bind-to core --report-bindings ../topomulticore_2D.exe < ../compiled/input8
done >> ./strong/8_MPI1cores2D.txt









export ACC_NUM_CORES=8

for j in {1..6}
do
    mpirun -np 1 --map-by ppr:1:node:pe=8 --bind-to core --report-bindings ../topoMULTI_alongX.exe < ../compiled/input8
done >> ./strong/1_MPI8coresX.txt

export ACC_NUM_CORES=4

for j in {1..6}
do
    mpirun -np 2 --map-by ppr:1:socket:pe=4 --bind-to core --report-bindings ../topoMULTI_alongX.exe < ../compiled/input8
done >> ./strong/2_MPI4coresX.txt

export ACC_NUM_CORES=2

for j in {1..6}
do
    mpirun -np 4 --map-by ppr:2:socket:pe=2 --bind-to core --report-bindings ../topoMULTI_alongX.exe < ../compiled/input8
done >> ./strong/4_MPI2coresX.txt

export ACC_NUM_CORES=1

for j in {1..6}
do
    mpirun -np 8 --map-by ppr:4:socket:pe=1 --bind-to core --report-bindings ../topoMULTI_alongX.exe < ../compiled/input8
done >> ./strong/8_MPI1coresX.txt









export ACC_NUM_CORES=8

for j in {1..6}
do
    mpirun -np 1 --map-by ppr:1:node:pe=8 --bind-to core --report-bindings ../topoMULTI_alongY.exe < ../compiled/input8
done >> ./strong/1_MPI8coresY.txt

export ACC_NUM_CORES=4

for j in {1..6}
do
    mpirun -np 2 --map-by ppr:1:socket:pe=4 --bind-to core --report-bindings ../topoMULTI_alongY.exe < ../compiled/input8
done >> ./strong/2_MPI4coresY.txt

export ACC_NUM_CORES=2

for j in {1..6}
do
    mpirun -np 4 --map-by ppr:2:socket:pe=2 --bind-to core --report-bindings ../topoMULTI_alongY.exe < ../compiled/input8
done >> ./strong/4_MPI2coresY.txt

export ACC_NUM_CORES=1

for j in {1..6}
do
    mpirun -np 8 --map-by ppr:4:socket:pe=1 --bind-to core --report-bindings ../topoMULTI_alongY.exe < ../compiled/input8
done >> ./strong/8_MPI1coresY.txt