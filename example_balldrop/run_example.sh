#!/bin/sh

echo `date`
echo "Run example of ball-drop impact."

mpirun -np 1 ../bin/swpc_3d.x -i ./input.inf
