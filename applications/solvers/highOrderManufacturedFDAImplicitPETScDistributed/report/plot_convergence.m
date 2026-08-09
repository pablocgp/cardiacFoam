function plot_convergence()
% Backwards-compatible wrapper. The convergence plots are now split into
%   plot_convergence_2D.m  (2D tutorial root only, fits over N=10..100 and
%                           N=70..100)
%   plot_convergence_3D.m  (three 3D tutorial roots, fits over N=10..30 and
%                           N=20..30)
% Both share the rendering routine in renderConvergence.m. Running either
% one alone is enough to refresh the figures for that dimension; this
% wrapper runs both for convenience.

plot_convergence_2D();
plot_convergence_3D();
end
