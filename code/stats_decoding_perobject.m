function stats_decoding_perobject()

% ----- run stats split by face/object/word: central condition -------
% all conditions versus chance (one-tailed v 0.5)
% comparisons across each hemisphere (two-tailed v 0)

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results and collate
fprintf('Loading data\n')

% get files
fileList = dir('results/sub-*_decoding_equal.mat');
nofiles = dir('results/sub-*_half_decoding.mat'); % don't include split half decoding files
filenames = {fileList.name};
idx = ismember(filenames,{nofiles.name});
files = filenames(~idx);
clust = {'L_elecs' 'R_elecs'};

% conds
objconds = {'objects','faces','words'};
stimnum.objects = [0:7 12:23];
stimnum.faces = 8:11;
stimnum.words = 24:2:46;

% collate
res_cell={};
cc = clock();mm='';
for f=1:length(files)
    fn = sprintf('results/%s',files{f});
    load(fn,'res');
    res = cosmo_slice(res,contains(res.sa.cond,'Central'));

    % set objtype and objnum as decoding between exemplars of each category
    % (e.g., word1 vs word2 = word category
    res.sa.objtype = cell(size(res.sa.clust));
    res.sa.objnum = zeros(size(res.sa.clust));
    for o = 1:length(objconds)
        idx = find(ismember(res.sa.cond1,stimnum.(objconds{o}))&...
            ismember(res.sa.cond2,stimnum.(objconds{o})));
        res.sa.objtype(idx) = objconds(o);
        res.sa.objnum(idx) = o;
    end

    % now make combos
    res.sa.objcombo1 = cell(size(res.sa.clust));
    res.sa.objcombo2 = cell(size(res.sa.clust));
    for o = 1:length(objconds)
        idx1 = find(ismember(res.sa.cond1,stimnum.(objconds{o})));
        idx2 = find(ismember(res.sa.cond2,stimnum.(objconds{o})));
        res.sa.objcombo1(idx1) = objconds(o);
        res.sa.objcombo2(idx2) = objconds(o);
    end

    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end
res_collate = cosmo_stack(res_cell);

res_collate.sa.combnum = zeros(size(res_collate.sa.objcombo1));
idx_f1 = contains(res_collate.sa.objcombo1,'face');
idx_f2 = contains(res_collate.sa.objcombo2,'face');
idx_o1 = contains(res_collate.sa.objcombo1,'object');
idx_o2 = contains(res_collate.sa.objcombo2,'object');
idx_w1 = contains(res_collate.sa.objcombo1,'word');
idx_w2 = contains(res_collate.sa.objcombo2,'word');

res_collate.sa.combnum(idx_f1&idx_f2) = 1;
res_collate.sa.combnum(idx_w1&idx_w2) = 2;
res_collate.sa.combnum(idx_o1&idx_o2) = 3;
res_collate.sa.combnum((idx_f1&idx_o2)|(idx_f2&idx_o1)) = 4;
res_collate.sa.combnum((idx_f1&idx_w2)|(idx_f2&idx_w1)) = 5;
res_collate.sa.combnum((idx_w1&idx_o2)|(idx_w2&idx_o1)) = 6;
combos = {'Face' 'Word' 'Object' 'FaceObject' 'FaceWord' 'WordObject'};

% average across all pairs
x=0;
reso=struct();
for p = 1:length(files)
    for cl = 1:2
        for cn = 1:6
            x=x+1;
            idx = find(res_collate.sa.subject==p&res_collate.sa.combnum==cn&contains(res_collate.sa.clust,clust{cl}));
            reso.samples(x,:) = mean(res_collate.samples(idx,:),1);
            reso.sa.subject(x,1) = p;
            reso.sa.combnum(x,1) = cn;
            reso.sa.clust{x,1} = clust{cl};
            fprintf('sub%d, cl%d, cn%d: %d rows\n', p, cl, cn, length(idx));
            nrows(p,cl,cn) = length(idx);
        end
    end
end

%% set up stats structure
clust = {'L_elecs' 'R_elecs'};

stats = struct();
stats.format = 'stats.objtype.electrodecluster';
stats.timevect = res_collate.a.fdim.values{1};

%% run stats
fprintf('Computing stats\n')

% collate data per object
for o = 1:length(combos)
    %% decoding per hemisphere
    for cl = 1:length(clust)

        % get means and do stats
        fprintf('Computing stats for %s, %s cluster\n',combos{o},clust{cl})

        dat = cosmo_slice(reso,strcmp(reso.sa.clust,clust{cl})&reso.sa.combnum==o);

        % get data details and stats for the group
        x = dat.samples;
        s = struct();
        s.eleccluster = clust{cl};
        s.objpair = combos(unique(dat.sa.combnum));
        s.n = size(x,1);
        s.mu = mean(x,1);
        s.mu_all = x;
        s.se = std(x)./sqrt(s.n);

        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

        % group mean onset & peak
        [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
        [~,peak] = max(s.mu); % peak
        ot(mo==0) = []; % remove "onsets" with max=0
        s.onset = stats.timevect(ot);
        s.peak = stats.timevect(peak);

        stats.(combos{o}).(clust{cl}) = s;
        fprintf('BFs done for %s, %s cluster\n',combos{o},clust{cl})

    end

    %% calculate hemispheric differences
    datdiff = stats.(combos{o}).(clust{2}).mu_all - stats.(combos{o}).(clust{1}).mu_all;

    % get data details and stats for the group difference
    s = struct();
    s.comparison = 'LH-RH difference';
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

    stats.(combos{o}).LHRHdiff = s;
end


fprintf('Saving\n')
save('results/stats_decoding_perobs.mat','stats');
fprintf('Done\n')


