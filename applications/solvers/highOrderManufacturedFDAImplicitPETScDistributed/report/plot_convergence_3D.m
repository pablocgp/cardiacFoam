function plot_convergence_3D()
% Plot 3D Picard/diagonalIion/JFNK convergence blocks for CN/consistent runs.

close all; clc;

resultRoot = '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitPETSc/results/example0/3D_run_20260522_235518';
renderPETScConvergenceBlocks('3D', resultRoot, [10 30], [20 30]);
end
