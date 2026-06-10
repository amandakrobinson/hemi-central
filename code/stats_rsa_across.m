function stats_rsa_across()

% assess similarity of representations in the left and right hemispheres
% use decoding performed on two halves of the data to make 2x RDMs
% correlate RDM from hemisphere 1 (half 1) to hemisphere 2 (half 2) and vice versa
% (do on separate halves to avoid inflation of correlation due to
% correlated decoding performance from noise)
% use spearman correlation
% take mean of two correlations LHh1:RHh2 and LHh2:RHh1

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results
fprintf('Loading data\n')
files = dir('results/sub-*_half_decoding.mat');
res_cell={};
cc = clock();mm='';
for f=1:length(files)
    fn = sprintf('results/%s',files(f).name);
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end
res_all = cosmo_stack(res_cell);

timevect = res_all.a.fdim.values{1};

%% computing left-right rsa hemisphere correlations
fprintf('Computing stats\n')

conds = {'Central'};
clust = {'L_elecs' 'R_elecs'};

stats = struct();
stats.dims = 'LHtime_RHtime_participant';
stats.measure = 'Spearman correlations';
stats.conds = conds;
stats.electrodeclusters = clust;

res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

across = zeros(length(timevect),length(timevect),length(files));
within = zeros(length(timevect),length(timevect),length(clust),length(files));
for p = 1:length(files)
    fprintf('Computing correlations for participant %d\n',p)

    res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

    %% across: correlate LH half 1 to RH half 2 and LH half 2 to RH half 1
    r = nan(length(timevect),length(timevect),length(clust));
    for h = 1:2 % for each split-half cross-validation

        % left hemisphere
        dat1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{1})&res_p.sa.half==h);

        % right hemisphere
        dat2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{2})&res_p.sa.half==(mod(h,2)+1));

        % get data details and stats for the group
        x1 = dat1.samples;
        x2 = dat2.samples;

        % get correlations
        r(:,:,h) = corr(x1,x2,'type','spearman'); % rows x cols is LH x RH

    end

    % take mean of h1 and h2
    across(:,:,p) = mean(r,3);


    %% correlate within hemis: LH half 1 to LH half 2, and RH half 1 to RH half 2
    r = nan(length(timevect),length(timevect),length(clust));
    for hemi = 1:2 % for each hemi

        % left hemisphere
        half1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{hemi})&res_p.sa.half==1);

        % right hemisphere
        half2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{hemi})&res_p.sa.half==2);

        % get data details and stats for the group
        x1 = half1.samples;
        x2 = half2.samples;

        % get correlations
        r(:,:,hemi) = corr(x1,x2,'type','spearman'); % rows x cols is LH x RH

    end

    % take mean of h1 and h2
    within(:,:,:,p) = r;


end
stats.across = across;
stats.within = within;
fprintf('Saving\n')
save('results/stats_across_hemi_correlations.mat','stats');

%% plot to check
figure(1);clf
subplot(1,3,1)
imagesc(timevect,timevect,mean(across,3),[-.1 .3])
axis square
title('LH-RH correlation')
subplot(1,3,2)
imagesc(timevect,timevect,mean(within(:,:,1,:),4),[-.1 .3])
axis square
title('LH split-half correlation')
subplot(1,3,3)
imagesc(timevect,timevect,mean(within(:,:,2,:),4),[-.1 .3])
axis square
title('RH split-half correlation')


%% now do mean correlation by left-right hemi offset
idxs = -100:100;
offsetcalcs = struct();

corr_tw = [];
for n = 1:size(across,3) % for each participant

    dn = across(:,:,n);

    % restrict to times 0-500ms
    idx = find(timevect>=0&timevect<=500);
    dn = dn(idx,idx);

    % now get mean correlations for off diagonals
    for j = 1:length(idxs)
        corr_tw(j,n) = mean(diag(dn,idxs(j)));
    end

end

mu = mean(corr_tw,2);
se = std(corr_tw,[],2)/sqrt(size(corr_tw,2));

mu_posneg = [];
for a = 1:100
    mu_posneg(a,1) = mu(idxs==a);
    mu_posneg(a,2) =  mu(idxs==(a*-1));
    xval(a,1) = a;
    xval(a,2) = a*-1;
end

% stats
bf = [];
for a = 1:100
    posdat = corr_tw(idxs==a,:);
    negdat = corr_tw(idxs==(a*-1),:);
    bf(a) = bayesfactor_R_wrapper(posdat-negdat,...
        'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
end

offsetcalcs.directional.offsettimes = idxs;
offsetcalcs.directional.correlation = corr_tw;
offsetcalcs.directional.mu = mu;
offsetcalcs.directional.se = se;

offsetcalcs.absolute.offsettimes = xval;
offsetcalcs.absolute.mu = mu_posneg;
offsetcalcs.absolute.bf = bf;

stats.offsetcalcs = offsetcalcs;

save('results/stats_across_hemi_correlations.mat','stats');

%% now do stats on correlations

rsa.Central = stats.data; % lh x rh

combinedstats=struct();
combinedstats.dims = 'leftTime_rightTime_participant';
combinedstats.measure = 'Spearman correlations';
combinedstats.electrodeclusters = clust;
fprintf('Pre-saving\n')
save('results/stats_across_hemi_correlations.mat','combinedstats','stats','timevect');

dat = rsa.Central;
clear across
across.n = size(dat,3);
across.tv = res_all.a.fdim.values{1};

for t = 1:size(dat,1) % for each LH hemi time point

    fprintf('Computing stats for timepoint %d of %d\n',t,size(dat,1))

    x = squeeze(dat(t,:,:))';

    across.mu(t,:) = mean(x);
    across.dat(t,:,:) = x';
    across.se(t,:) = std(x)./sqrt(across.n);

    % calculate bayesfactors
    across.bf(t,:) = bayesfactor_R_wrapper(x',...
        'returnindex',2,'verbose',false,...
        'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

end

combinedstats.data = across;

timevect = across.tv;
fprintf('Saving\n')
save('results/stats_across_hemi_correlations.mat','combinedstats','stats','timevect');
fprintf('Done\n')


