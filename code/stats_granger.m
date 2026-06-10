function stats_granger()

% function to do granger causality across hemispheres
% partial correlations between LH(t) and RH(t+15) while partialling out RH(t),
% showing what information came from the other hemi
% use eg RDM from odd sequences correlated with RDM from even

% correlate using spearman correlation

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox (Personal)/MATLAB/CoSMoMVPA/mvpa')
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

%% computing granger
fprintf('Computing stats\n')

clust = {'L_elecs' 'R_elecs'};

stats = struct();

res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

for p = 1:length(files)
    fprintf('Computing correlations for participant %d\n',p)

    res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

    %% correlate across hemispheres with time delay, partialling out other hemi
    RtoL = NaN(length(timevect),7,2);
    LtoR = NaN(length(timevect),7,2);
    for ha = 1:2 % for each half of data -do twice

        % left hemi, one half of sequences
        dat1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{1})&res_p.sa.half==ha);

        % right hemi, other half of sequences
        dat2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,clust{2})&res_p.sa.half==(mod(ha,2)+1));

        % get data details and stats for the group
        x1 = dat1.samples;
        x2 = dat2.samples;

        % correlate hemispheres per time point
        lags = 1:4:30;
        for l = 1:length(lags)
            lag = lags(l);
            for t = 31:length(timevect) % from -80 ms
                RtoL(t,l,ha) = partialcorr(x1(:,t),x2(:,t-lag),x1(:,t-lag),'type','Spearman');
                LtoR(t,l,ha) = partialcorr(x2(:,t),x1(:,t-lag),x2(:,t-lag),'type','Spearman');
            end
        end

    end
    stats.granger.RtoL(:,:,p) = mean(RtoL,3);
    stats.granger.LtoR(:,:,p) = mean(LtoR,3);

end

fprintf('Saving\n')
save('results/stats_granger.mat','stats','timevect');
