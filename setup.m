function cfg = setup()
%SETUP Paths and every setting of the phase-transition experiment.

root = fileparts(mfilename('fullpath'));
lib = fullfile(root, 'lib');

%% Paths
addpath(fullfile(lib, 'RPCA'), fullfile(lib, 'PROPACK'));
% PREEEDM and RGNMR (optional, see README). 
base = fullfile(lib, 'Baselines');
if isfolder(fullfile(base, 'PREEEDM'))
    addpath(fullfile(base, 'PREEEDM', 'PREEEDM', 'PREEEDM'));
    addpath(fullfile(lib, 'PREEEDM_patch'));
end
if isfolder(fullfile(base, 'RGNMR'))
    addpath(fullfile(base, 'RGNMR', 'MATLAB'), fullfile(base, 'RGNMR', 'MATLAB', 'utils', 'RGNMR_utils'));
    addpath(fullfile(lib, 'RGNMR_patch'));
end
addpath(fullfile(lib, 'core'));
if ~(exist('haltonset', 'file') && exist('pdist', 'file') && exist('randsample', 'file'))
    error('setup:Toolbox', 'The Statistics and Machine Learning Toolbox is required (haltonset, pdist, randsample).');
end
if ~exist('wthresh', 'file')        % Wavelet Toolbox missing: same hard thresholding in plain MATLAB
    addpath(fullfile(lib, 'core', 'compat'));
end
if isfolder(fullfile(pwd, 'PROPACK'))
    warning('setup:Folder', ['Run from the main folder: get_mu_kappa (lib/RPCA) switches to a different SVD ' ...
        'routine when the current folder has a PROPACK subfolder.']);
end

%% Data
cfg.dims = [2 3];                  % dimensions d
cfg.m_values = 10:5:70;            % numbers of anchors
cfg.alpha_values = 0:0.05:0.5;     % outlier fractions
cfg.n_trials = 50;                 % trials per (d, m, alpha)
cfg.T = 500;                       % total number of points
cfg.box = 100;
cfg.outlier_scale = 1;
cfg.base_seed = 22281015;

%% Methods
cfg.methods = {'SREDG'};           % the paper also runs 'PREEEDM' and 'RGNMR' (clone them first, see README)
cfg.sredg = struct('mu_factor', 1.1, 'tol', 1e-14, 'gamma', 0.9, 'max_iter', 500, 'trimming', false);
cfg.preeedm = struct('max_iter', 2000);   % PREEEDM's internal cap
% RGNMR: (m+n) x (m+n) matrix with E and F observed, the true number of corrupted entries given.
cfg.rgnmr = struct('init_option', 0, 'max_outer_iter', 50, 'max_inner_iter', 2000, 'inner_init_tol', 1e-15, ...
    'stop_relRes', 1e-15, 'stop_relDiff', -1, 'LSQR_smart_tol', 0, 'verbose', 1);

%% Figure
cfg.fig.thresholds = [1];     % threshold for recovery
cfg.fig.colormap = 'gray';         % 0 = black, 1 = white
cfg.fig.nodata_color = [0.957 0.839 0.839];
cfg.fig.panel_size = [260 230];    % pixels per panel (width, height)
cfg.fig.font_size = 9;
cfg.fig.formats = {'png'};         
cfg.fig.resolution = 300;          % dpi of the png

%% Run
cfg.workers = 0;                   % 0 = no parallelization, or the number of workers (e.g. 10)
cfg.results_dir = fullfile(root, 'results');   % each run makes results/run_<date>_<time>/
end
