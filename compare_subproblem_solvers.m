% Compare the Lanczos and conjugate gradient subproblem solvers of ARC
% running on a given problem in Manopt.
%
% First version: May 2, 2019
%
% Bryan Zhu, Nicolas Boumal
% https://github.com/bwzhu97/arc_subproblem

clear; clc;

% Fix randomness once and for all
rng(2018);

cd example_problems;
addpath(genpath(pwd()));

%% Run ARC on a problem
problem = rotation_synchronization(3, 50, .75);
init = problem.M.rand();
solver = struct('solver', @arc_sub_comparison, 'theta', 0.05, 'tolgradnorm', 1e-8);
solver.statsfun = statsfunhelper(statscounters({'hesscalls', 'gradhesscalls'}));
[x, cost, info] = manoptsolve(problem, init, solver);
substats_l = {info.lanczos};
substats_cg = {info.cg};
sigmas = [info.sigma];
iterations = length(sigmas);

cd ..;

idstring = datestr(now(), 'mmm_dd_yyyy_HHMMSS');

%% Plot results
cd outputs;
addpath(genpath(pwd()));

subplot_rows = 3;
subplot_cols = 3;
enough_points = 6;
plots = subplot_rows * subplot_cols;

outputs = {'gradnorms', 'func_values'};
yscale = {'log', 'linear'};
axisnames = {'Gradient norm', 'Model cost'};
noutputs = numel(outputs);

for output = 1 : noutputs
    figure(output);
    clf;
    set(gcf, 'Color', 'w');
    output_name = outputs{output};
    
    I = 1;
    filled_plots = 0;
    while I < iterations && filled_plots < plots
        output_l = substats_l{I+1}.(output_name);
        output_cg = substats_cg{I+1}.(output_name);
        if length(output_l) < enough_points && length(output_cg) < enough_points
            I = I + 1;
            continue;
        end

        subplot(subplot_rows, subplot_cols, filled_plots+1);
        plot_title = sprintf('Iteration %d, sigma = %.2e', I, sigmas(I));
        title(plot_title);

        hold all;
        plot(output_cg, ...
             'DisplayName', 'ARC-CG', ...
             'Marker', '.', 'MarkerSize', 15);
        plot(output_l, ...
             'DisplayName', 'ARC-L', ...
             'Marker', '.', 'MarkerSize', 15);
        hold off;

        set(gca, 'XScale', 'linear');
        set(gca, 'YScale', yscale{output});
        if filled_plots == 0
            legend('show');
        end
        if mod(filled_plots, subplot_cols) == 0
            ylabel(axisnames{output});
        end
        if floor(filled_plots / subplot_cols) == subplot_rows - 1
            xlabel('Inner iteration #');
        end
        grid on;

        I = I + 1;
        filled_plots = filled_plots + 1;
    end
    figname = sprintf('compare_lanczos_cg_%s_%s', idstring, output_name);
    savefig([figname, '.fig']);
    %pdf_print_code(gcf, [figname, '.pdf'], 14);
end


cd ..;