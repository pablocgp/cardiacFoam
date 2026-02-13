#!/bin/bash

cd "$1"
source $HOME/OpenFOAM-v2312/etc/bashrc

echo "🧹 Cleaning case"
./Allclean

echo "🚀 Running parallel"
./Allrun

echo "✍️ Touching case.foam for paraview"
touch case.foam


