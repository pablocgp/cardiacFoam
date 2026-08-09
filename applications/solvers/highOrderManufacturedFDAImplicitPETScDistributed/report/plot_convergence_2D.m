function plot_convergence_2D()
% Plot 2D Picard/diagonalIion/JFNK convergence blocks for CN/consistent runs.

close all; clc;

resultRoot = '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitPETSc/results/example0/2D_run_20260522_034403';
renderPETScConvergenceBlocks('2D', resultRoot, [10 100], [70 100]);
end
