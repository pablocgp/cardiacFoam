function plot_convergence_2D()
% Plot 2D cell-centred L2 and reconstructed G2 convergence figures for the
% implicit manufactured FDA JfNK/Picard/diagonalIion runs.
%
% One figure file is generated for each nonlinear method, time step, time
% scheme, mass matrix, variable, and error metric. The legend reports the
% fitted convergence orders over N=10..100 (general) and N=70..100
% (asymptotic), together with the accumulated runtime and peak-RSS cost over
% the mesh sweep represented by each curve.
%
% Companion script plot_convergence_3D.m generates the 3D figures from the
% three separate tutorial roots used for the 3D sweep.
%
% NOTE: result directories still use the legacy token "hoStates". In the
% current solver this token denotes the high-order quadrature level used for
% Iion and for the local ODE state values stored at the Iion integration
% points.

close all; clc;
dimTag = '2D';
nRangeAll = [10 100];
nRangeFine = [70 100];

methodSources = struct( ...
    'Picard', '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK/results/example0/2D', ...
    'diagonalIion', '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK/results/example0/2D', ...
    'JFNK', '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK/results/example0/2D' ...
);

renderConvergence(dimTag, methodSources, nRangeAll, nRangeFine);
end
