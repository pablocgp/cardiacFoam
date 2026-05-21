function plot_convergence_3D()
% Plot 3D cell-centred L2 and reconstructed G2 convergence figures for the
% implicit manufactured FDA JfNK/Picard/diagonalIion runs.
%
% The 3D sweep is split across three tutorial roots (one per nonlinear
% method); the "_1" sibling directories may contain in-progress runs and
% are searched as fallback sources. The legend reports the fitted
% convergence orders over the general N=10..30 mesh range and over the
% asymptotic N=20..30 range, plus the accumulated runtime and peak-RSS cost
% of the mesh sweep represented by each curve.
%
% Companion script plot_convergence_2D.m generates the 2D figures from the
% single 2D tutorial root.
%
% NOTE: result directories still use the legacy token "hoStates". In the
% current solver this token denotes the high-order quadrature level used for
% Iion and for the local ODE state values stored at the Iion integration
% points.

close all; clc;
dimTag = '3D';
nRangeAll = [10 30];
nRangeFine = [20 30];

methodSources = struct( ...
    'Picard', {{ ...
        '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK_Picard/results/example0/3D', ...
        '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK_Picard_1/results/example0/3D' ...
    }}, ...
    'diagonalIion', {{ ...
        '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK_diagonalIion/results/example0/3D', ...
        '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK_diagonalIion_1/results/example0/3D' ...
    }}, ...
    'JFNK', {{ ...
        '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK/results/example0/3D' ...
    }} ...
);

renderConvergence(dimTag, methodSources, nRangeAll, nRangeFine);
end
