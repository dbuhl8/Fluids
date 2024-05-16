#/bin/bash

mpirun --use-hwthread-cpus -np 20 shearSolve >> output.txt
