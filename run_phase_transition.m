function run_dir = run_phase_transition(cfg)
%RUN_PHASE_TRANSITION Recovery of SREDG, PREEEDM and RGNMR over a grid of anchors m and outlier fractions alpha.
%
%   run_phase_transition();          % settings from setup.m
%   cfg = setup(); cfg.methods = {'SREDG'}; run_phase_transition(cfg);
%
% For every dimension d, number of anchors m and outlier fraction alpha in cfg, the run draws cfg.n_trials random
% instances, runs every method in cfg.methods on each, and records the RMSE of the recovered points.

if nargin < 1, cfg = setup(); end
% PREEEDM and RGNMR need extra setup (see README)
probe = {'PREEEDM', 'dcroot'; 'RGNMR', 'check_early_convergence'};
for i = 1:size(probe, 1)
    if ismember(probe{i, 1}, cfg.methods) && ~exist(probe{i, 2}, 'file')
        error('run_phase_transition:Missing', 'Clone %s into lib/Baselines/%s first (see README).', probe{i, 1}, probe{i, 1});
    end
end

%% Results folder: settings and console log
run_dir = fullfile(cfg.results_dir, ['run_' char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'))]);
mkdir(run_dir);
fid = fopen(fullfile(run_dir, 'settings.json'), 'w');
fprintf(fid, '%s\n', jsonencode(cfg, 'PrettyPrint', true));
fclose(fid);
diary(fullfile(run_dir, 'log.txt'));
stop_diary = onCleanup(@() diary('off'));

%% Parallel pool (cfg.workers > 0) and one computational thread per trial, as in the paper runs
nw = 0;
if cfg.workers > 0
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= cfg.workers
        delete(pool);
        cluster = parcluster('Processes');
        cluster.NumWorkers = max(cluster.NumWorkers, cfg.workers);  
        pool = parpool(cluster, cfg.workers);
    end
    nw = pool.NumWorkers;
end
old_threads = maxNumCompThreads(1);
restore_threads = onCleanup(@() maxNumCompThreads(old_threads));

%% Main loop: one (d, m, alpha) cell at a time; the trials of a cell run in parallel
n_cells = numel(cfg.dims) * numel(cfg.m_values) * numel(cfg.alpha_values);
fprintf('%s | %s | d = %s | %d cells x %d trials\n', run_dir, strjoin(cfg.methods, ', '), ...
    mat2str(cfg.dims), n_cells, cfg.n_trials);
file = fullfile(run_dir, 'per_trial.csv');
results = table;
k = 0;
t0 = tic;
for d = cfg.dims
    for m = cfg.m_values
        for alpha = cfg.alpha_values
            k = k + 1;

            % Run the trials of this cell
            rows = cell(cfg.n_trials, 1);
            parfor (t = 1:cfg.n_trials, nw)
                rows{t} = run_trial(cfg, d, m, alpha, t);
            end

            % Save them right away, so a stopped run keeps its finished cells
            R = struct2table(vertcat(rows{:}));
            writetable(R, file, 'WriteMode', 'append');
            results = [results; R]; %#ok<AGROW>

            % Progress: recovery rate (RMSE < 1) of each method in this cell, and time left
            rate = cellfun(@(mm) sprintf('%s %3.0f%%', mm, 100 * mean(R.rmse(strcmp(R.method, mm)) < 1)), ...
                cfg.methods, 'UniformOutput', false);
            fprintf('[%3d/%d] d = %d, m = %2d, alpha = %.2f | RMSE < 1: %s | about %s left\n', k, n_cells, d, m, ...
                alpha, strjoin(rate, ', '), char(duration(0, 0, (n_cells - k) * toc(t0) / k, 'Format', 'hh:mm:ss')));
        end
    end
end

%% Summary table and one phase diagram per recovery threshold
writetable(summarize(results, cfg.fig.thresholds), fullfile(run_dir, 'summary.csv'));
for thr = cfg.fig.thresholds
    out = fullfile(run_dir, ['phase_rmse_lt_' strrep(sprintf('%g', thr), '.', 'p')]);
    plot_phase_transition(results, thr, out, cfg.fig);
end
fprintf('Done in %s: %d outcomes, %d failed. Results in %s\n', ...
    char(duration(0, 0, toc(t0), 'Format', 'hh:mm:ss')), height(results), nnz(~strcmp(results.status, 'ok')), run_dir);
end


%% ============================================================================================================
function rows = run_trial(cfg, d, m, alpha, trial)
% One instance and every method on it. Returns one row per method.

% The seed depends only on (d, m, alpha, trial), so any sub-grid repeats exactly the same instances
seed = mod(cfg.base_seed + 1000003 * (d - 2) + 7919 * m + 104729 * round(1000 * alpha) + trial, 2^31 - 1);
inst = make_instance(cfg, d, m, alpha, seed);

% Every method starts from the random-generator state right after the instance, whatever the method order
state = rng;
for j = 1:numel(cfg.methods)
    rng(state);
    rows(j, 1) = run_method(cfg.methods{j}, inst, cfg, trial); %#ok<AGROW>
end
end


function inst = make_instance(cfg, d, m, alpha, seed)
% One instance: T points in [-box, box]^d, the first m of them anchors.
% E = anchor-anchor squared distances (clean), F = anchor-target squared distances (with outliers).
rng(seed, 'twister');

% Points: anchors from a scrambled Halton set, targets uniform at random
n = cfg.T - m;
hs = scramble(haltonset(d, 'Skip', 1e3, 'Leap', 1e2), 'RR2');
P1 = -cfg.box + 2 * cfg.box * net(hs, m);           % m x d anchors
P2 = -cfg.box + 2 * cfg.box .* rand(d, n);          % d x n targets
P = [P1', P2];

% Squared distances
D = squareform(pdist(P')).^2;
E = D(1:m, 1:m);
F = D(1:m, m+1:end);

% Outliers: a fraction alpha of the entries of F, uniform on [-s*mean|F|, s*mean|F|] with s = outlier_scale
supp = randsample(m * n, round(alpha * m * n), false);
S_range = cfg.outlier_scale * mean(abs(F(:)));
S_temp = 2 * S_range * rand(m, n) - S_range;
S = zeros(m, n);
S(supp) = S_temp(supp);
Fobs = F + S;

% AccAltProj parameters for SREDG; the incoherence mu is estimated from the observed F
r = d + 2;                                          % rank of the squared-distance matrix
mu = cfg.sredg.mu_factor * get_mu_kappa(Fobs, r);
para = struct('mu', mu, 'beta_init', r * sqrt(mu(1) * mu(end)) / sqrt(m * n), ...
    'beta', r * sqrt(mu(1) * mu(end)) / (4 * sqrt(m * n)), 'trimming', cfg.sredg.trimming, ...
    'tol', cfg.sredg.tol, 'gamma', cfg.sredg.gamma, 'max_iter', cfg.sredg.max_iter);

inst = struct('d', d, 'm', m, 'alpha', alpha, 'seed', seed, 'P', P, 'E', E, 'Fobs', Fobs, ...
    'n_outliers', nnz(S), 'rank', r, 'para', para);
end


function row = run_method(method, inst, cfg, trial)
% One method on one instance: estimated points Y (T x d), then the RMSE against the true points.
d = inst.d; m = inst.m; r = inst.rank; E = inst.E; Fobs = inst.Fobs;
row = struct('d', d, 'm', m, 'alpha', inst.alpha, 'trial', trial, 'seed', inst.seed, 'method', method, ...
    'rmse', NaN, 'iterations', NaN, 'converged', false, 'hit_cap', false, 'seconds', NaN, ...
    'status', 'ok', 'error', '');
timer = tic;
try
    switch method
        case 'SREDG'
            % Clean F with AccAltProj, complete the Gram matrix (Nystrom), take its top-d eigenvectors.
            % The iteration count comes from AccAltProj's printed log.
            para = inst.para;
            [txt, X] = quietly(@() SREDG(E, Fobs, r, @AccAltProj, para));
            tok = regexp(txt, 'Total (\d+) iteration', 'tokens', 'once');
            row.converged = ~isempty(tok);
            if row.converged, row.iterations = str2double(tok{1}); else, row.iterations = para.max_iter; end
            row.hit_cap = row.iterations >= para.max_iter;
            [V, L] = eigs(X, d, 'lm');
            Y = V * sqrt(L);

        case 'PREEEDM'
            rng(inst.seed + 9000000, 'twister');           % PREEEDM's own seed for this trial
            n = size(Fobs, 2);

            % Distance matrix (not squared); negative observed squared distances are clipped to 0
            Dr = zeros(m + n);
            Dr(1:m, 1:m) = sqrt(E);
            Dr(1:m, m+1:end) = sqrt(max(Fobs, 0));
            Dr(m+1:end, 1:m) = Dr(1:m, m+1:end)';

            % Observed pairs: anchor-anchor and anchor-target, zero distances included
            W = false(m + n);
            W(1:m, 1:m) = ~eye(m);
            W(1:m, m+1:end) = true;
            W(m+1:end, 1:m) = true;

            % Anchor positions: classical MDS of E
            J = eye(m) - ones(m) / m;
            B = -J * E * J / 2; B = (B + B') / 2;
            [V, L] = eig(B, 'vector');
            [L, idx] = sort(L, 'descend');
            A = (V(:, idx(1:d)) * diag(sqrt(max(L(1:d), 0))))';

            % Solve; the targets are o.X 
            pars = struct('m', m, 'PP', [A zeros(d, n)], 'draw', 0, 'observation_mask', W);
            [~, o] = quietly(@() PREEEDM_zeroaware(Dr, d, pars));
            Y = [A o.X]';
            row.iterations = o.iterations;
            row.hit_cap = o.iterations >= cfg.preeedm.max_iter;
            row.converged = ~row.hit_cap && o.ErrEig < .01 && o.ErrObj < log(nnz(W)) * 1e-4;

        case 'RGNMR'
            % Full (m+n) x (m+n) squared-distance matrix with E and F observed, target-target missing
            n = size(Fobs, 2); N = m + n;
            X = zeros(N);
            X(1:m, 1:m) = E;
            X(1:m, m+1:N) = Fobs;
            X(m+1:N, 1:m) = Fobs';
            mask = false(N);
            mask(1:m, :) = true;
            mask(:, 1:m) = true;
            [ii, jj] = find(mask);
            omega = [ii, jj];
            k = 2 * inst.n_outliers;                        

            % Solve
            opts = cfg.rgnmr;
            rng(inst.seed + 9000000, 'twister');           % RGNMR's own seed for this trial
            [txt, Dhat, flag] = quietly(@() RGNMR(X, omega, r, k, opts));
            tok = regexp(txt, '\[INSIDE GNMR\] iter\s+(\d+)', 'tokens');
            if ~isempty(tok), row.iterations = str2double(tok{end}{1}); end
            row.converged = logical(flag);
            row.hit_cap = row.iterations >= opts.max_outer_iter;
            Y = edm_to_points(Dhat, d);

        otherwise
            error('Unknown method %s', method);
    end

    % RMSE after centering and orthogonal Procrustes alignment
    assert(isreal(Y) && all(isfinite(Y(:))), 'Nonreal or nonfinite coordinates');
    [~, row.rmse] = quietly(@() Compute_RMSE(inst.P', Y));
catch err
    row.status = 'failed';
    row.error = err.message;
end
row.seconds = toc(timer);
end


function [txt, varargout] = quietly(f) %#ok<INUSD> (f is called inside evalc)
% Calls f() without printing; txt is the text it would have printed.
[txt, varargout{1:nargout-1}] = evalc('f()');
end


function Q = edm_to_points(Dhat, d)
% Points from a completed squared-distance matrix: symmetrize, zero the diagonal, classical MDS.
N = size(Dhat, 1);
D = (Dhat + Dhat') / 2;
D(1:N+1:end) = 0;
J = eye(N) - ones(N) / N;
B = -.5 * J * D * J; B = (B + B') / 2;
[V, L] = eig(B, 'vector');
[L, idx] = sort(L, 'descend');
assert(all(L(1:d) > 0), 'Fewer than d positive Gram eigenvalues');
Q = V(:, idx(1:d)) * diag(sqrt(L(1:d)));
end


function S = summarize(T, thresholds)
% One row per (method, d, m, alpha)
[g, keys] = findgroups(T(:, {'method', 'd', 'm', 'alpha'}));
S = keys;
S.trials = splitapply(@numel, T.rmse, g);
for thr = thresholds
    S.(['recovered_rmse_lt_' strrep(sprintf('%g', thr), '.', 'p')]) = splitapply(@(x) mean(x < thr), T.rmse, g);
end
S.median_rmse = splitapply(@(x) median(x, 'omitnan'), T.rmse, g);
S.failed = splitapply(@(s) nnz(~strcmp(s, 'ok')), T.status, g);
S.hit_cap = splitapply(@nnz, T.hit_cap, g);
S.median_seconds = splitapply(@median, T.seconds, g);
end
