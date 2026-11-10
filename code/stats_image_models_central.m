function stats_image_models_central()

% run RSA: correlate visual statistics RDMs with neural RDMs from left and
% right hemispheres (separately) while partialling out the other hemisphere

doonsets=1;

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
end

%% load image models
load('results/RSAmodels.mat','m_cielab','m_lum','m_rect','m_sflow','m_sfhigh','m_edge');
mask = ones(size(m_cielab));
mask = tril(mask,-1);

mods = [m_edge(mask==1) m_cielab(mask==1) m_lum(mask==1) ...
    m_sflow(mask==1) m_sfhigh(mask==1) m_rect(mask==1)];
modelnames = {'edge' 'colour' 'luminance' 'lowSF' 'highSF'...
    'rectilinearity'};

%% get neural rdms

load('results/rdms.mat','alldat')
timevect = alldat.time;
clust = {'L_elecs' 'R_elecs'};

%% correlate neural data with models
fprintf('Computing correlation analyses\n')

clear s
for p = 1:20
    fprintf('Computing RSA for participant %d\n',p)

    %% partial correlations between each hemisphere and image models, controlling for other hemi

    % left hemisphere
    x1 = alldat.(clust{1})(:,:,p);

    % right hemisphere
    x2 = alldat.(clust{2})(:,:,p);

    % straight correlations
    s.corr_L(:,:,p) = corr(x1,mods,'type','Spearman');
    s.corr_r(:,:,p) = corr(x2,mods,'type','Spearman');

    % group means for noise ceilings
    low_lh = mean(alldat.L_elecs(:,:,setdiff(1:20,p)),3); % mean of group except p
    low_rh = mean(alldat.R_elecs(:,:,setdiff(1:20,p)),3);
    high_lh = mean(alldat.L_elecs,3); % mean of group including p
    high_rh = mean(alldat.R_elecs,3);

    % partial correlations
    for t = 1:length(timevect)
        % participant partial correlation - one hemi with image models, partialling out other hemi
        s.partialcorr_L(t,:,p) = partialcorr(x1(:,t),mods,x2(:,t),'type','Spearman');
        s.partialcorr_r(t,:,p) = partialcorr(x2(:,t),mods,x1(:,t),'type','Spearman');

        % lower noise ceiling
        s.partialceil_lowerL(t,:,p) = partialcorr(x1(:,t),low_lh(:,t),x2(:,t),'type','Spearman');
        s.partialceil_lowerR(t,:,p) = partialcorr(x2(:,t),low_rh(:,t),x1(:,t),'type','Spearman');
        
        % upper noise ceiling
        s.partialceil_upperL(t,:,p) = partialcorr(x1(:,t),high_lh(:,t),x2(:,t),'type','Spearman');
        s.partialceil_upperR(t,:,p) = partialcorr(x2(:,t),high_rh(:,t),x1(:,t),'type','Spearman');
    end
end
%%
dat = s;
fprintf('Saving\n')
save('results/stats_imagecorrs.mat','dat','timevect','mods','modelnames', '-v7.3');

%% now do stats on correlations

stats = struct();
stats.timevect = timevect;
stats.clust = clust;

contrast = {'partialcorr_L', 'partialcorr_r'};

combs = nchoosek(1:20,20-2); %leave-2-participants-out for jackknifing

for c = 1:length(contrast)
    for m = 1:length(modelnames)
        fprintf('Computing stats for %s, %s\n',contrast{c}, modelnames{m})
        s=struct();
        x = squeeze(dat.(contrast{c})(:,m,:));

        mu = mean(x,2);
        se = std(x,[],2)./sqrt(size(x,2));
        s.mu_all = x;
        s.mu = mu;
        s.se = se;

        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(x,...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

        % group mean onset & peak
        [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
        [~,peak] = max(s.mu); % peak
        ot(mo==0) = []; % remove "onsets" with max=0
        s.onset = stats.timevect(ot);
        s.peak = stats.timevect(peak);

        %% calculate onsets and peaks by jackknifing
        tv = timevect<400;
        if doonsets
            fprintf('Calculating onsets for %s\n',contrast{c})

            % use jackknife method and calculate bayesfactors to get onsets
            x=s.mu_all';
            xa = [];
            for i=1:size(combs,1)
                xa = [xa x(combs(i,:),tv)];
            end
            tic
            b = bayesfactor_R_wrapper(xa',...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
            s.bf_jackknife = reshape(b,[],size(combs,1));
            fprintf('jackknife BFs done for %s\n',contrast{c})
            toc

            % get onsets and peaks per jackknife
            onset_jack = [];peak_jack=[];maxon=[];
            for b=1:size(combs,1)
                [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
                [~,peak_jack(b)] = max(mean(x(combs(b,:),:),1));
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
        end

        stats.(contrast{c}).(modelnames{m}) = s;
        fprintf('Onsets done\n')

        fprintf('Saving\n')
        save('results/stats_imagecorrs.mat','stats','dat','timevect','mods','modelnames', '-v7.3');
        fprintf('Done\n')
    end
end

%% now do stats on hemisphere differences

contrast = {'partialcorr'};

for c = 1:length(contrast)
    for m = 1:length(modelnames)
        fprintf('Computing stats for %s, %s\n',contrast{c}, modelnames{m})

        RH = squeeze(dat.(sprintf('%s_r',contrast{c}))(:,m,:));
        LH = squeeze(dat.(sprintf('%s_L',contrast{c}))(:,m,:));

        x=RH-LH;
        mu = mean(x,2);
        se = std(x,[],2)./sqrt(size(x,2));
        stats.(sprintf('diff_%s',contrast{c})).(modelnames{m}).mu = mu;
        stats.(sprintf('diff_%s',contrast{c})).(modelnames{m}).se = se;

        % calculate bayesfactors
        stats.(sprintf('diff_%s',contrast{c})).(modelnames{m}).bf = bayesfactor_R_wrapper(x,...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

    end
end

fprintf('Saving\n')
save('results/stats_imagecorrs.mat','stats','dat','timevect','mods','modelnames', '-v7.3');
fprintf('Done\n')

