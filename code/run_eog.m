%run eog analysis

% Eye-movement proxy analysis (reviewer point 2), central condition.
% No dedicated EOG was recorded, so frontal electrodes are used as a proxy:
%   HEOG_F  = F8 - F7          (positive = rightward gaze)
% For each trial: baseline-corrected (-200 to 0 ms) mean amplitude and
% variance over 0-400 ms. Summaries over all trials and per stimulus
% (stimC 0:23 and 24:2:48)
%
% Uses the non-CSD data: the Laplacian suppresses the far-field ocular
% signals this analysis relies on.

%% for each participant

for subjectnr=1:20

    fn    = sprintf('../hemifield-eeg/data/derivatives/eeglab/sub-%02i_task-rsvp_continuous.set',subjectnr);
    outfn = sprintf('results/sub-%02i_eog.mat',subjectnr);

    % load
    fprintf('loading %s\n',fn);
    EEG = pop_loadset(fn);
    eventlist = readtable(sprintf('../hemifield-eeg/data/sub-%02i/eeg/sub-%02i_task-rsvp_events.csv',subjectnr,subjectnr));

    % create individual trial epochs
    EEG = pop_epoch(EEG, {'E  1'}, [-.2 .4]);
    assert(EEG.trials == height(eventlist), 'sub-%02i: %i epochs, %i events', ...
        subjectnr, EEG.trials, height(eventlist));

    % central condition, selected stimuli only
    idx = find(eventlist.condition==1);
    events = eventlist(idx,:);
    EEG = pop_select(EEG,'trial',idx);
    EEG = pop_rmbase(EEG,[-200 0]);

    timevect = EEG.times;
    vartime = find(EEG.times>=0);

    % compute
    chanidx = [find(ismember({EEG.chanlocs(:).labels},'F7')) find(ismember({EEG.chanlocs(:).labels},'F8'))];
    HEOG = squeeze(EEG.data(chanidx(2),:,:)-EEG.data(chanidx(1),:,:)); % time x trial

    % fix stim number for words
    its = 25:2:47; % italics indices
    cits=ismember(events.stimC,its);
    events.stimC(cits) = events.stimC(cits)-1;

    stim = unique(events.stimC);
    x=0;
    for s = 1:numel(stim)
        x=x+1;
        idx = events.stimC==stim(s);
        erp(x,:,subjectnr) = mean(HEOG(:,idx), 2);                 % stim x time x subj
        v(x,subjectnr)     = mean(var(HEOG(vartime,idx), 0, 1));   % stim x subj    end
    end
end

save('results/heog.mat','erp','v','timevect')

%% objects, faces, words
% stimnum.objects = [0:7 12:23];
% stimnum.faces = 8:11;
% stimnum.words = 24:2:46;

objs = erp([1:8 13:24],:,:); % stim x time x ps
faces = erp(9:12,:,:);
words = erp(25:36,:,:);
v_objs = v([1:8 13:24],:); % stim x ps
v_faces = v(9:12,:);
v_words = v(25:36,:);

bfv=bayesfactor_R_wrapper(mean(v_objs,1)-mean(v_words,1),...
'returnindex',1,'verbose',false,...
'args','mu=0,rscale="medium"');   % BF10, paired
fprintf('Difference between trial variance of object and word trials BF=%.03g\n',bfv)
bfv=bayesfactor_R_wrapper(mean(v_objs,1)-mean(v_faces,1),...
'returnindex',1,'verbose',false,...
'args','mu=0,rscale="medium"');   % BF10, paired
fprintf('Difference between trial variance of object and face trials BF=%.03g\n',bfv)
bfv=bayesfactor_R_wrapper(mean(v_words,1)-mean(v_faces,1),...
'returnindex',1,'verbose',false,...
'args','mu=0,rscale="medium"');   % BF10, paired
fprintf('Difference between trial variance of word and face trials BF=%.03g\n',bfv)

cols = tab10(3);


figure(1);clf;
set(gcf,'Position',[37 332 840 840])
subplot(2,1,1);hold on
fill([0 100 100 0],[-2 -2 2 2],[.7 .7 .7],'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
text(8,1.7,'Stimulus on','Fontsize',18)

meanobjs = mean(objs,1);
mu = mean(meanobjs,3);
se = std(meanobjs,[],3)/sqrt(size(meanobjs,3));
fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],cols(1,:),'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
plot(timevect,mu,'Color',cols(1,:),'LineWidth',2)

meanf = mean(faces,1);
mu = mean(meanf,3);
se = std(meanf,[],3)/sqrt(size(meanf,3));
fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],cols(3,:),'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
plot(timevect,mu,'Color',cols(3,:),'LineWidth',2)

meanwords = mean(words,1);
mu = mean(meanwords,3);
se = std(meanwords,[],3)/sqrt(size(meanwords,3));
fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],cols(2,:),'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
plot(timevect,mu,'Color',cols(2,:),'LineWidth',2)


title('HEOG by category')
legend({'Objects' 'Faces' 'Words'})
set(gca,'FontSize',20)
ylim([-2 2])
ylabel('HEOG proxy, F8-F7 difference (\muV)')
xlabel('Time (ms)')



% calculate and plot BFs
for d = 1:3
    if d == 1
        owdiff = squeeze(meanobjs-meanf);
        name = 'Objects versus faces';
    elseif d == 2
        owdiff = squeeze(meanobjs-meanwords);
        name = 'Objects versus words';
    else
        owdiff = squeeze(meanf-meanwords);
        name = 'Faces versus words';
    end

    bf(:,d)=bayesfactor_R_wrapper(owdiff,...
        'returnindex',1,'verbose',false,...
        'args','mu=0,rscale="medium"');   % BF10, paired
    fprintf('Max BF for %s is %.03g\n',name,max(bf(:,d)))

    sp=subplot(6,1,3+d);hold on;
    plot(timevect,1+0*timevect,'k-');
    idx = [bf(:,d)<1/10,1/10<bf(:,d) & bf(:,d)<10,bf(:,d)>10]';
    col = [0 0 0];
    cob = [.5 .5 .5;1 1 1;col];
    for i=1:3
        x = timevect(idx(i,:));
        y = bf(idx(i,:),d);
        if ~isempty(x)
            stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
        end
    end
    sp.YScale='log';
    ylim([10.^-2 10.^2])
    xlim([-200 400])
    sp.YTick = 10.^([-2 0 2]);
    sp.YTickLabel = {'0.01' '1' '100'};
    sp.TickDir = 'out';
    sp.XTickLabel = {};
    set(sp,'FontSize',22)
    sp.Position(4) = sp.Position(4)-.01;
    text(-170,30,name,'FontSize',18)
    
    if d == 2
        ylabel('Bayes Factors')
    end
end
xlabel('Time (ms)')

fn = 'figures/FigureS2_HEOG_3cats';
set(gcf, 'Renderer', 'painters')
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

names = {'Objects versus faces' 'Objects versus words' 'Faces versus words'};
for a = 1:3
    less1 = sum(bf(:,a)<1)/length(timevect);
    less3 = sum(bf(:,a)<3)/length(timevect);
    fprintf('%s: BF<1 %.2f, BF<3 %.2f\n',names{a},less1*100,less3*100)
end
