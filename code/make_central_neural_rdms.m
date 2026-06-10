function make_central_neural_rdms()

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
end
load('results/imageorder.mat')

%% stack results
fprintf('Loading data\n')
ps = 1:20;
res_cell={};
cc = clock();mm='';
for f=1:length(ps)
    fn = sprintf('results/sub-%02i_decoding.mat',ps(f));
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(ps),sprintf('%i/%i',f,length(ps)),mm);
end
res_all = cosmo_stack(res_cell);

%% collate left-right RDMs

clust = {'L_elecs' 'R_elecs'};

alldat = struct();
alldat.imagecombs = combnk(0:1:35,2);
alldat.time = res_all.a.fdim.values{1};


res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

for p = 1:length(ps)
    fprintf('Computing RDMs for participant %d\n',p)

    res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

    % left hemisphere
    dat1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{1}));

    % right hemisphere
    dat2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{2}));

    % get data details and stats for the group
    alldat.(clust{1})(:,:,p) = dat1.samples;
    alldat.(clust{2})(:,:,p) = dat2.samples;

end
fprintf('Saving\n')
save('results/rdms.mat','alldat');


alldat.ims = [dat1.sa.cond1 dat2.sa.cond2];
alldat.imnames = [dat1.sa.im1 dat2.sa.im2];
alldat.measure = 'pairwise_decoding_accuracy';

fprintf('Saving\n')
save('results/rdms.mat','alldat');


