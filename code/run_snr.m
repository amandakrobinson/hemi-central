function run_snr(subjectnr)

% Signal quality (SNR) analysis for the left and right electrode clusters,
% central condition, following Luck et al. (2021) Psychophysiology 58:e13793
% "Standardized measurement error (SME)".
%
% For each participant and cluster (LH, RH) computes:
%   - the ERP (mean over all trials and the 6 cluster electrodes)
%   - time-window mean amplitude for "P1" and "N1" windows (pos and neg
%   components after stim onset)
%   - analytic SME for those mean amplitudes (SD of single-trial scores / sqrt(n))
%
% Uses the same preprocessed data (CSD/Laplacian) and electrode clusters as
% run_decoding_equal.m, so signal quality is assessed on the data that were decoded.

rng(subjectnr);

%% paths
if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
end

%% settings
fn    = sprintf('../hemifield-eeg/data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa_CSD.mat',subjectnr);
outfn = sprintf('results/sub-%02i_snr.mat',subjectnr);

elecsel = [14:16 44:46; 18:20 48:50];   % same clusters as run_decoding.m
clustnames = {'LH' 'RH'};

windows = struct( ...
    'name',  {'P1',      'N1'}, ...
    'range', {[80 120],  [160 200]});   % ms

%% load
fprintf('loading %s\n',fn);
load(fn,'ds')

% central condition, valid central stimuli only
ds = cosmo_slice(ds, ds.sa.condition==1 & ds.sa.stimC>-1);

% trials x channels x time
[arr, dimlabels, dimvalues] = cosmo_unflatten(ds);
assert(strcmp(dimlabels{1},'chan') && strcmp(dimlabels{2},'time'), 'unexpected feature dimensions');
timevect = dimvalues{2};
timevect = timevect(:)';

%% compute
rows = {};
erp = struct();
for c = 1:2
    % average over the 6 cluster electrodes -> trials x time
    x = squeeze(mean(arr(:, elecsel(c,:), :), 2));
    n  = size(x,1);

    erp.(clustnames{c}) = mean(x,1);

    for w = 1:numel(windows)
        win = timevect>=windows(w).range(1) & timevect<=windows(w).range(2);
        wt  = timevect(win);

        % --- time-window mean amplitude and analytic SME ---
        trialscores = mean(x(:,win), 2);
        meanamp = mean(trialscores);
        aSME    = std(trialscores) / sqrt(n);
        snr     = abs(meanamp) / aSME;


        rows(end+1,:) = {subjectnr, clustnames{c}, windows(w).name, ...
            n, meanamp, aSME, snr};

    end
end

snr = cell2table(rows, 'VariableNames', {'subject','cluster','window', ...
    'ntrials','meanamp','aSME','SNR'});

if ~isfolder('results'), mkdir('results'); end
save(outfn, 'snr', 'erp', 'timevect', 'windows')
fprintf('saved %s\n', outfn);

end