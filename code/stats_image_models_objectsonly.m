function stats_image_models_objectsonly()

% run RSA: correlate visual statistics RDMs with neural RDMs from left and right hemispheres
% objects only
% partial correlations: one hemi with visual stats, accounting for other
% hemi

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

%% sort out subsets of decoding pairs
objconds = {'objects','words'};
stimnum.objects = [0:7 12:23];
stimnum.words = 24:2:46;

objidx = zeros(length(mods),1);
combs = nchoosek([0:23 24:2:46],2); % objects, faces and words (36)
for o = 1:2
    idx = ismember(combs(:,1),stimnum.(objconds{o}))&...
        ismember(combs(:,2),stimnum.(objconds{o}));
    objidx(idx) = o;
end

clust = {'L_elecs' 'R_elecs'};

%% get neural rdms

load('results/rdms.mat','alldat')
timevect = alldat.time;

%% correlate neural data with models
fprintf('Computing correlation analyses\n')

clear s
for p = 1:size(alldat.L_elecs,3)

    for o = 1:length(objconds)
        fprintf('Computing RSA for participant %d\n',p)

    %% partial correlations LH ~ visual models | RH (and reverse)

        % left hemisphere
        x1 = alldat.(clust{1})(:,:,p);

        % right hemisphere
        x2 = alldat.(clust{2})(:,:,p);

        % straight correlations
        s.corr_L(:,:,o,p) = corr(x1(objidx==o,:),mods(objidx==o,:),'type','Spearman');
        s.corr_r(:,:,o,p) = corr(x2(objidx==o,:),mods(objidx==o,:),'type','Spearman');

        % partial correlations
        for t = 1:length(timevect)
            s.partialcorr_L(t,:,o,p) = partialcorr(x1(objidx==o,t),mods(objidx==o,:),x2(objidx==o,t),'type','Spearman');
            s.partialcorr_r(t,:,o,p) = partialcorr(x2(objidx==o,t),mods(objidx==o,:),x1(objidx==o,t),'type','Spearman');
        end
    end
end
dat = s;
fprintf('Saving\n')
save('results/stats_imagecorrs_perobjs.mat','dat','timevect','mods','modelnames', '-v7.3');

%% now do stats on correlations

stats = struct();
stats.timevect = timevect;
stats.clust = clust;

contrast = fieldnames(dat);

for c = 1:length(contrast)
    for m = 1:length(modelnames)
        for o = 1:length(objconds)
            fprintf('Computing stats for %s, %s, %s\n',contrast{c}, modelnames{m},objconds{o})

            x = squeeze(dat.(contrast{c})(:,m,o,:));

            mu = mean(x,2);
            se = std(x,[],2)./sqrt(size(x,2));
            stats.(contrast{c}).(objconds{o}).(modelnames{m}).mu = mu;
            stats.(contrast{c}).(objconds{o}).(modelnames{m}).se = se;
            stats.(contrast{c}).(objconds{o}).(modelnames{m}).mu_all = x;

            % calculate bayesfactors
            stats.(contrast{c}).(objconds{o}).(modelnames{m}).bf = bayesfactor_R_wrapper(x,...
                'returnindex',2,'verbose',false,...
                'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
        end
    end
end

fprintf('Saving\n')
save('results/stats_imagecorrs_perobjs.mat','stats','dat','timevect','mods','modelnames','objconds', '-v7.3');
fprintf('Done\n')


%% now do stats on hemisphere differences

contrast = {'corr' 'partialcorr'};

for c = 1:length(contrast)
    for m = 1:length(modelnames)

        for o = 1:length(objconds)
            fprintf('Computing stats for %s, %s, %s\n',contrast{c}, modelnames{m},objconds{o})

            RH = squeeze(dat.(sprintf('%s_r',contrast{c}))(:,m,o,:));
            LH = squeeze(dat.(sprintf('%s_L',contrast{c}))(:,m,o,:));

            x=RH-LH;
            mu = mean(x,2);
            se = std(x,[],2)./sqrt(size(x,2));
            stats.(sprintf('diff_%s',contrast{c})).(objconds{o}).(modelnames{m}).mu = mu;
            stats.(sprintf('diff_%s',contrast{c})).(objconds{o}).(modelnames{m}).se = se;

            % calculate bayesfactors
            stats.(sprintf('diff_%s',contrast{c})).(objconds{o}).(modelnames{m}).bf = bayesfactor_R_wrapper(x,...
                'returnindex',2,'verbose',false,...
                'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
        end
    end
end

fprintf('Saving\n')
save('results/stats_imagecorrs_perobjs.mat','stats','dat','timevect','mods','modelnames','objconds', '-v7.3');
fprintf('Done\n')

%%
for o = 1:length(objconds)

    omods = mods(objidx==o,:);

    fig = figure(o);clf;
    fig.Position=[1 1 1400 800];

    for m = 1:size(mods,2)
        subplot(3,3,m)
        imagesc(squareform(omods(:,m)))
        axis square
        title(sprintf('%s, %s only', modelnames{m},objconds{o}))
    end
end

