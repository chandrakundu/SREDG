% RUN_DEMO  Quick phase-transition run


cfg = setup();

%% Quick setup
cfg.T = 500;                       % total number of points
cfg.dims = [2];                    % dimensions d
cfg.m_values = 50:5:60;            % numbers of anchors
cfg.alpha_values = 0.1:0.05:0.2;   % outlier fractions
cfg.n_trials = 5;                  % number of trials
cfg.methods = {'SREDG','PREEEDM','RGNMR'};      % methods to run, to compare with PREEEDM, RGNMR; follow instructions in the readme.
cfg.workers = 0;                   % 0 = no parallelization, or the number of workers (e.g. 10)

run_dir = run_phase_transition(cfg);
S = readtable(fullfile(run_dir, 'summary.csv'));
disp(S(:, {'method', 'd', 'm', 'alpha', 'trials', 'recovered_rmse_lt_1', 'failed'}));
