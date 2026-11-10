function stats_decoding_central()

% ----- run stats on decoding performance: central condition -------
% each hemisphere versus chance (one-tailed v 0.5)

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results
fprintf('Loading data\n')

% get files
fileList = dir('results/sub-*_decoding_equal.mat');
nofiles = dir('results/sub-*_half_decoding.mat'); % don't include split half decoding files
filenames = {fileList.name};
idx = ismember(filenames,{nofiles.name});
files = filenames(~idx);

doonsets=1;
combs = nchoosek(1:length(files),(length(files)-2)); %leave-2-participants-out for jackknifing

%% collate
res_cell={};
cc = clock();mm='';
for f=1:length(files)
    fn = sprintf('%s/%s',fileList(1).folder,files{f});
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end
res_all = cosmo_stack(res_cell);

% average across all pairs
res_all = cosmo_average_samples(res_all,'split_by',{'clust','cond','subject','vf'});

%% set up stats
fprintf('Computing stats\n')

conds = {'Central'};
clust = {'L_elecs' 'R_elecs'};

stats = struct();
stats.format = 'stats.electrodecluster';
stats.timevect = res_all.a.fdim.values{1};

%% run

res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

for cl = 1:length(clust)

    %% get means and do stats
    fprintf('Computing stats for central, %s cluster\n',clust{cl})

    dat = cosmo_slice(res_cond,strcmp(res_cond.sa.clust,clust{cl}));

    % get data details and stats for the group
    x = dat.samples;
    s = struct();
    s.condition = conds{1};
    s.eleccluster = clust{cl};
    s.n = size(x,1);
    s.mu = mean(x,1);
    s.mu_all = x;
    s.se = std(x)./sqrt(s.n);

    % calculate bayesfactors
    s.bf = bayesfactor_R_wrapper(x',...
        'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

    s.tv = res_all.a.fdim.values{1};

    % group mean onset & peak
    [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
    [~,peak] = max(s.mu); % peak
    ot(mo==0) = []; % remove "onsets" with max=0
    s.onset = stats.timevect(ot);
    s.peak = stats.timevect(peak);

    stats.(clust{cl}) = s;
    fprintf('BFs done for %s cluster \n',clust{cl})

    %% calculate onsets and peaks by jackknifing
    if doonsets == 1
        fprintf('Calculating onsets for %s cluster \n',clust{cl})

        % use jackknife method and calculate bayesfactors to get onsets
        x = [];
        for i=1:size(combs,1)
            x = [x s.mu_all(combs(i,:),:)];
        end
        tic
        b = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');
        s.bf_jackknife = reshape(b,[],size(combs,1));
        fprintf('jackknife BFs done for %s cluster \n',clust{cl})
        toc

        % get onsets and peaks per jackknife
        onset_jack = [];peak_jack=[];maxon=[];
        for b=1:size(s.bf_jackknife,2)
            [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
            [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:),1));
        end
        onset_jack(maxon==0) = []; % remove "onsets" with max=0

        % convert onsets and peaks to time (rather than tp)
        s.onsetjk = stats.timevect(onset_jack);
        if isempty(onset_jack)
            s.onsetci=[];
        else
            ci = prctile(onset_jack,[2.5 97.5]);
            s.onsetci = stats.timevect([floor(ci(1)) ceil(ci(2))]);
        end
        s.peakjk = stats.timevect(peak_jack);
        s.peakci = stats.timevect(round(prctile(peak_jack,[2.5 97.5])));

        stats.(clust{cl}) = s;
        fprintf('Onsets done\n')
    end
end

%% RH-LH hemispheric differences

datdiff = stats.R_elecs.mu_all - stats.L_elecs.mu_all; % RH-LH

s = struct();
s.comparison = 'RH-LH difference';
s.n = size(datdiff,1);
s.mu = mean(datdiff,1);
s.mu_all = datdiff;
s.se = std(datdiff)./sqrt(s.n);

% calculate bayesfactors
s.bf = bayesfactor_R_wrapper(datdiff',...
    'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

% group mean onset & peak
[mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
[~,peak] = max(s.mu); % peak
ot(mo==0) = []; % remove "onsets" with max=0
s.onset = stats.timevect(ot);
s.peak = stats.timevect(peak);

if doonsets == 1
    fprintf('Calculating onsets for hemi diff\n')

    % use jackknife method and calculate bayesfactors to get onsets
    x = [];
    for i=1:size(combs,1)
        x = [x s.mu_all(combs(i,:),:)];
    end
    tic
    b = bayesfactor_R_wrapper(x',...
        'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
    s.bf_jackknife = reshape(b,[],size(combs,1));
    fprintf('jackknife BFs done for %s cluster\n',clust{cl})
    toc

    % get onsets and peaks per jackknife
    onset_jack = [];peak_jack=[];maxon=[];
    for b=1:size(s.bf_jackknife,2)
        [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
        [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:),1));
    end
    onset_jack(maxon==0) = []; % remove "onsets" with max=0

    % convert onsets and peaks to time (rather than tp)
    s.onsetjk = stats.timevect(onset_jack);
    if isempty(onset_jack)
        s.onsetci=[];
    else
        conf = prctile(onset_jack,[2.5 97.5]);
        s.onsetci = stats.timevect([floor(conf(1)) ceil(conf(2))]);
    end
    s.peakjk = stats.timevect(peak_jack);
    s.peakci = stats.timevect(round(prctile(peak_jack,[2.5 97.5])));

    fprintf('Onsets done\n')
end

stats.RHLHdiff = s;

%%
fprintf('Saving\n')
save('results/stats_decoding_central.mat','stats');
fprintf('Done\n')


