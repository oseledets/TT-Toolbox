function [y,swp]=mvk3(a,x,eps,varargin)
%  Approximates the tt matvec W*X in the norm via
%   the DMRG iterations.
%   default varargin are the following:
%   kickrank = 5;
%   dropsweeps = 1; % garbage - for quasi-wedderburn
%   % dropsweeps2 = 10; % garbage
%   % d_pow_trunc = 1.5; % garbage
%   ddpow = 0.1; % stepsize for d-power in truncations
%   ddrank = 1; % stepsize for additional rank
%   d_pow_check = 0; % d-power for checking the convergence
%   bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
%   top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank
%   verb = 1; % 0 - silent, 1 - sweep information, 2 - block information
%   rmax=1000; % maximal rank
%   nswp=25; % maximal number of sweep
% wrapper for the tt_mvk3!

[y,swp]=tt_mvk3(core(a),core(x),eps,varargin{:});
y=tt_tensor(y);