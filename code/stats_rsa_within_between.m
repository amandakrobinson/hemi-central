function stats_rsa_within_between()

% function to correlate data halves for a given hemisphere
% eg RDM from odd sequences correlated with RDM from even
% measure of data consistency for contra vs ipsi

% then compare to correlations across hemispheres

% correlate using spearman correlation

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
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

%% computing within hemi reliability: half-half rsa correlations within hemispheres/VF conditions
fprintf('Computing stats\n')

clust = {'L_elecs' 'R_elecs'};

stats = struct();

res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

for p = 1:length(files)
    fprintf('Computing correlations for participant %d\n',p)

    res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

    %% correlate within a hemisphere: consistency
    for hemi = 1:2 % for each hemisphere

        % first half (odd/even)
        dat1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{hemi})&res_p.sa.half==1);

        % second half (odd/even)
        dat2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{hemi})&res_p.sa.half==2);

        % get data details and stats for the group
        x1 = dat1.samples;
        x2 = dat2.samples;

        % correlate halves per pairs of time points
        r=corr(x1,x2,'type','Spearman');

        stats.(clust{hemi}).mu_all(:,:,p) = r; % half1 x half2
    end

    %% correlate across hemispheres: shared information
    acr = zeros(length(timevect),length(timevect),2);
    for ha = 1:2 % for each half of data -do twice

        % left hemi, one half of sequences
        dat1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{1})&res_p.sa.half==ha);

        % right hemi, other half of sequences
        dat2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{2})&res_p.sa.half==(mod(ha,2)+1));

        % get data details and stats for the group
        x1 = dat1.samples;
        x2 = dat2.samples;

        % correlate hemispheres per time point
        acr(:,:,ha) = corr(x1,x2,'type','Spearman');

    end
    stats.across.raw(:,:,:,p) = acr; % lh x rh x half x participant
    stats.across.mu_all(:,:,p) = mean(acr,3);

end

fprintf('Saving\n')
save('results/stats_withinacross.mat','stats','timevect');

%% do stats for correlations different from zero

con1 = [clust 'across'];

for x = 1:length(con1) % for each hemisphere/across

    dat = stats.(con1{x}).mu_all;
    clear s
    s.mu_all = dat;
    s.n = size(dat,3);
    s.mu = mean(dat,3);
    s.se = std(dat,[],3)./sqrt(s.n);

    % calculate bayesfactors
    s.bf = NaN(length(timevect),length(timevect));
    for t = 1:size(s.mu,1)
        fprintf('Contrast %s: calculating stats for timepoint %d\n',con1{x},t)
        s.bf(t,:) = bayesfactor_R_wrapper(squeeze(dat(t,:,:)),...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
    end
    stats.(con1{x}) = s;
    save('results/stats_withinacross.mat','stats','timevect');

end

%% now do differences between left-across and right-across
for ci = 1:length(clust)
    dat1 = stats.(clust{ci}).mu_all; % contra/ipsi
    dat2 = stats.across.mu_all; % across
    clear s
    dat = dat1-dat2;
    s.mu_all = dat;
    s.n = size(dat,3);
    s.mu = mean(dat,3);
    s.se = std(dat,[],3)./sqrt(s.n);

    % calculate bayesfactors
    s.bf = NaN(length(timevect),length(timevect));
    for t = 1:size(s.mu,1)
        fprintf('%s within vs across: calculating stats for timepoint %d\n',con1{ci},t)

        s.bf(t,:) = bayesfactor_R_wrapper(squeeze(dat(t,:,:)),...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
    end
    stats.withinVacross.(clust{ci}) = s;
end


fprintf('Saving\n')
save('results/stats_withinacross.mat','stats','timevect');
fprintf('Done\n')

%% now do left-right differences

dat1 = stats.L_elecs.mu_all; % lh
dat2 = stats.R_elecs.mu_all; % rh
clear s
dat = dat2-dat1;
s.mu_all = dat;
s.n = size(dat,3);
s.mu = mean(dat,3);
s.se = std(dat,[],3)./sqrt(s.n);

% calculate bayesfactors
s.bf = NaN(length(timevect),length(timevect));
for t = 1:size(s.mu,1)
    fprintf('RH minus LH: calculating stats for timepoint %d\n',t)

    s.bf(t,:) = bayesfactor_R_wrapper(squeeze(dat(t,:,:)),...
        'returnindex',2,'verbose',false,...
        'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
end
stats.LHRHdiff = s;


fprintf('Saving\n')
save('results/stats_withinacross.mat','stats','timevect');
fprintf('Done\n')

