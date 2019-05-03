% Compare the Lanczos and conjugate gradient subproblem solvers of ARC on
% a cubic subproblem model with a randomly generated matrix and vector in
% Manopt.
%
% First version: May 2, 2019
%
% Bryan Zhu, Nicolas Boumal
% https://github.com/bwzhu97/arc_subproblem

clear; clc;

% Fix randomness once and for all
rng(2018);

%% Run the subproblem solvers on the cubic model
dim = 500;
[problem, b] = generate_hard_case(dim);
x = problem.M.rand();
bnorm = problem.M.norm(x, b);
sigma = 0.001;

options = struct('theta', 0, 'verbosity', 0, 'maxiter_lanczos', 250, 'maxiter_cg', 250);
storedb = StoreDB(1);
key = storedb.getNewKey();

[~, ~, ~, ~, substats_l] = arc_lanczos(problem, x, b, bnorm, sigma, options, storedb, key);
[~, ~, ~, ~, substats_cg] = arc_conjugate_gradient(problem, x, b, bnorm, sigma, options, storedb, key);

idstring = datestr(now(), 'mmm_dd_yyyy_HHMMSS');

%% Plot results
cd outputs;
addpath(genpath(pwd()));
outputs = {'gradnorms', 'func_values'};
yscale = {'log', 'log'};
axisnames = {'Gradient norm', 'Model cost'};
noutputs = numel(outputs);

for output = 1 : noutputs
    figure(output);
    clf;
    set(gcf, 'Color', 'w');
    output_name = outputs{output};

    output_l = substats_l.(output_name);
    output_cg = substats_cg.(output_name);

    title('Hard case Lanczos vs CG');
    hold all;
    plot(output_cg, 'DisplayName', 'ARC-CG', 'Marker', '.', 'MarkerSize', 15);
    plot(output_l, 'DisplayName', 'ARC-L', 'Marker', '.', 'MarkerSize', 15);
    hold off;

    set(gca, 'XScale', 'linear');
    set(gca, 'YScale', yscale{output});
    legend('show');
    ylabel(axisnames{output});
    xlabel('Inner iteration #');
    grid on;

    figname = sprintf('hardcase_%s_%s', idstring, output_name);
    savefig([figname, '.fig']);
end
cd ..;

function [problem, b] = generate_hard_case(dim)
    Q = RandOrthMat(dim);
    eigenvals = sort(randn(dim, 1), 'ascend');
    lambda = diag(eigenvals);
    A = Q.' * lambda * Q;
    q1 = Q(1,:);
    rvec = randn(dim, 1);
    b = rvec - ((q1 * rvec) / norm(q1)) * q1.';
        
    problem.M = euclideanfactory(dim);
    problem.hess = @hess;
    function [h, store] = hess(~, xdot, store)
        h = A*xdot;
    end
end


