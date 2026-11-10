function make_central_neural_rdms()

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
end
% load('results/imageorder.mat')

%% stack results
fprintf('Loading data\n')
ps = 1:20;
res_cell={};
cc = clock();mm='';
for f=1:length(ps)
    fn = sprintf('results/sub-%02i_decoding_equal.mat',ps(f));
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(ps),sprintf('%i/%i',f,length(ps)),mm);
end
res_all = cosmo_stack(res_cell);

%% collate left-right RDMs

clust = {'L_elecs' 'R_elecs'};

alldat = struct();
% alldat.imagecombs = combnk(0:1:35,2);
alldat.time = res_all.a.fdim.values{1};

res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

for p = 1:length(ps)
    res_p = cosmo_slice(res_cond, res_cond.sa.subject == ps(p));
    dat1 = cosmo_slice(res_p, strcmp(res_p.sa.clust, clust{1}));
    dat2 = cosmo_slice(res_p, strcmp(res_p.sa.clust, clust{2}));

    pairs1 = [dat1.sa.cond1 dat1.sa.cond2];
    pairs2 = [dat2.sa.cond1 dat2.sa.cond2];
    if p == 1, pairs_ref = pairs1; end
    assert(isequal(pairs1, pairs2, pairs_ref), 'pair order mismatch for sub-%02i', ps(p))

    alldat.(clust{1})(:,:,p) = dat1.samples;
    alldat.(clust{2})(:,:,p) = dat2.samples;
end

% stimulus codes and their RDM positions (1-36 = stim001-036)
codes = [0:23 24:2:46];
[~, pos] = ismember(pairs_ref, codes);
assert(all(pos(:) > 0))
alldat.ims        = pairs_ref;      % actual stimulus codes
alldat.imagecombs = pos - 1;        % 0-35 positions, as before
alldat.imnames    = [dat1.sa.im1 dat1.sa.im2];
alldat.measure    = 'pairwise_decoding_accuracy';

save('results/rdms.mat', 'alldat');

