function stats_image_models_partial()

% run RSA: correlate visual statistics RDMs with neural RDMs from left and right hemispheres
% regression on full model (one hemi ~ other hemi + 6 image models)
% noise ceiling: regression on full model ~ group mean
% unique variance calculate by regression of one hemi ~ all but one model

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
res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,'Central'));

timevect = res_all.a.fdim.values{1};
clust = {'L_elecs' 'R_elecs'};

%% correlate neural data with models
fprintf('Computing regression analyses\n')
modnums = 1:length(modelnames);
clear s
for p = 1:20
    fprintf('Computing RSA for participant %d\n',p)
    res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

    %% partial out one hemisphere from another, then regress with image models

    % left hemi first half (odd/even)
    dat1 = cosmo_slice(res_p,strcmp(res_p.sa.clust,'L_elecs')&res_p.sa.half==1);

    % right hemi second half (odd/even)
    dat2 = cosmo_slice(res_p,strcmp(res_p.sa.clust,'R_elecs')&res_p.sa.half==2);


    % left hemisphere
    x1 = dat1.samples;

    % right hemisphere
    x2 = dat2.samples;

    % regress and get R2 vals
    for t = 1:length(timevect)

        % full model
        mdl_L = fitlm([mods,x2(:,t)],x1(:,t));
        mdl_R = fitlm([mods,x1(:,t)],x2(:,t));
        s.stats_all(t,:,p) = [mdl_L.Rsquared.Adjusted mdl_R.Rsquared.Adjusted];

        % without other hemi
        mdl_L_image = fitlm(mods,x1(:,t));
        mdl_R_image = fitlm(mods,x2(:,t));
        s.stats_image(t,:,p) = [mdl_L_image.Rsquared.Adjusted mdl_R_image.Rsquared.Adjusted];

        % without individual image models
        for m = 1:length(modelnames)
            mdl_Lmim = fitlm([mods(:,modnums(modnums~=m)),x2(:,t)],x1(:,t));
            mdl_Rmim = fitlm([mods(:,modnums(modnums~=m)),x1(:,t)],x2(:,t));
            s.stats_minusim(t,:,m,p) = [mdl_Lmim.Rsquared.Adjusted mdl_Rmim.Rsquared.Adjusted];
        end
        
    end

end

dat = s;
fprintf('Saving\n')
save('results/stats_regress.mat','dat','timevect','mods','modelnames', '-v7.3');

%% plot

subplot(2,3,1)
plot(timevect,mean(dat.stats_all,3))
subplot(2,3,2)
plot(timevect,mean(dat.stats_image,3))
subplot(2,3,3)
plot(timevect,mean(dat.stats_all-dat.stats_image,3))
title('unique other hemisphere')

% get unique each image stat
for p = 1:20
    for m = 1:6
        unmod(:,:,m,p) = dat.stats_all(:,:,p)-dat.stats_minusim(:,:,m,p);
    end
end
cols = tab10(7);
for h = 1:2
    subplot(2,3,3+h)
    hold on
    for m = 1:6
        d = squeeze(unmod(:,h,m,:));
        mu = mean(d,2);
        se = std(d,[],2)/sqrt(size(d,2));
        fill([timevect fliplr(timevect)],[mu'+se' fliplr(mu'-se')],cols(m,:),...
            'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off')
        plot(timevect,movmean(mu,10),'Color',cols(m,:),'LineWidth',2)
        title(clust{h})
        legend(modelnames)
        ylim([0 .06])
    end
end


figure
cols = viridis(3);
for m = 1:6
    subplot(2,3,m)
    hold on
    for h = 1:2
        d = squeeze(unmod(:,h,m,:));
        mu = mean(d,2);
        se = std(d,[],2)/sqrt(size(d,2));
        fill([timevect fliplr(timevect)],[mu'+se' fliplr(mu'-se')],cols(h,:),...
            'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off')
        plot(timevect,mu,'Color',cols(h,:),'LineWidth',2)
        title(modelnames{m})
        legend(clust)
        ylim([0 .06])
    end
end

%% get unique each image stat without other hemi
for p = 1:20
    for m = 1:6
        unmod(:,:,m,p) = dat.stats_image(:,:,p)-dat.stats_minusim(:,:,m,p);
    end
end

figure
cols = tab10(7);
for h = 1:2
    subplot(2,3,3+h)
    hold on
    for m = 1:6
        d = squeeze(unmod(:,h,m,:));
        mu = mean(d,2);
        se = std(d,[],2)/sqrt(size(d,2));
        fill([timevect fliplr(timevect)],[mu'+se' fliplr(mu'-se')],cols(m,:),...
            'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off')
        plot(timevect,movmean(mu,10),'Color',cols(m,:),'LineWidth',2)
        title(clust{h})
        legend(modelnames)
        ylim([0 .06])
    end
end


%% now do stats on correlations

stats = struct();
stats.timevect = timevect;
stats.clust = clust;

contrast = fieldnames(dat);
contrast = contrast(1:4); 

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
        if doonsets
            fprintf('Calculating onsets for %s\n',contrast{c})

            % use jackknife method and calculate bayesfactors to get onsets
            x=s.mu_all';
            xa = [];
            for i=1:size(combs,1)
                xa = [xa x(combs(i,:),:)];
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
                ci = prctile(onset_jack,[5 95]);
                s.onsetci = stats.timevect([floor(ci(1)) ceil(ci(2))]);
            end
            s.peakjk = stats.timevect(peak_jack);
            s.peakci = stats.timevect(round(prctile(peak_jack,[5 95])));
        end

        stats.(contrast{c}).(modelnames{m}) = s;
        fprintf('Onsets done\n')

        fprintf('Saving\n')
        save('results/stats_imagecorrs.mat','stats','dat','timevect','mods','modelnames', '-v7.3');
        fprintf('Done\n')
    end
end

%% now do stats on hemisphere differences

contrast = {'corr' 'partialcorr'};

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


%% neural-feature correlations, partialling out other models
clear s

for p = 1:20
    fprintf('Computing RSA for participant %d\n',p)

    %% correlation neural with one model, partialling out others

    % left hemisphere
    x1 = alldat.(clust{1})(:,:,p);

    % right hemisphere
    x2 = alldat.(clust{2})(:,:,p);

    % partial correlations
    for m = 1:length(modelnames)
        for t = 1:length(timevect)
            s.unique_L(t,m,p) = partialcorr(x1(:,t),mods(:,m),mods(:,setdiff(1:6,m)),'type','Spearman');
            s.unique_R(t,m,p) = partialcorr(x2(:,t),mods(:,m),mods(:,setdiff(1:6,m)),'type','Spearman');
        end
    end
end
partialdat = s;
fprintf('Saving\n')
save('results/stats_imagecorrs.mat','stats','dat','partialdat','timevect','mods','modelnames', '-v7.3');
