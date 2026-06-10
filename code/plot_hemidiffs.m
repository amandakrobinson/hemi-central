
tws = [90 120; 135 165; 210 240];

cols = tab10(6);

figure(1);clf;set(gcf,'Position',[1687 140 1457 900])

%% all stimuli
ylims = [-.15 .22];

load('/Users/uqarobi7/Dropbox/PROJECTS/hemiraptor/hemi-central/results/stats_imagecorrs.mat')
subplot(3,2,1)
plot(timevect,timevect*0,'LineWidth',0.5,'Color','k','HandleVisibility','off')
hold on
for m = 1:length(modelnames)
    diff.(modelnames{m}).mu = stats.partialcorr_r.(modelnames{m}).mu-stats.partialcorr_L.(modelnames{m}).mu;
    diff.(modelnames{m}).mu_all = stats.partialcorr_r.(modelnames{m}).mu_all-stats.partialcorr_L.(modelnames{m}).mu_all;

    plot(timevect,movmean(diff.(modelnames{m}).mu,20),'LineWidth',2,'Color',cols(m,:))
end
for t = 1:size(tws,1)
    fill([tws(t,:) fliplr(tws(t,:))],[-.15 -.15 .25 .25],[.7 .7 .7],'FaceAlpha',.2,'EdgeColor','none')
end
legend(modelnames)

for t = 1:3
    tw = tws(t,:);
    tw = find(timevect>=tw(1)&timevect<=tw(2));
    for m = 1:length(modelnames)
        mdat = mean(diff.(modelnames{m}).mu_all(tw,:),1);
        mu(t,m) = mean(mdat);
        se(t,m) = std(mdat)/sqrt(length(mdat));
    end
    timenames{t} = sprintf('%d to %d ms',tws(t,1),tws(t,2));
end
xlabel('Time (ms)')
xlim([-100 800])
set(gca,'FontSize',20)
ylabel('RH-LH difference')
title('All stimuli')
ylim(ylims)

subplot(3,2,2); hold on
nGroups = size(tws,1);
nBars = 6;
groupWidth = 0.7;
barWidth = groupWidth / nBars;
for i = 1:nBars
    % calculate x positions for this bar across both groups
    xPos = (1:nGroups) - groupWidth/2 + (i-0.5) * barWidth;
    
    % draw bars
    bar(xPos, mu(:,i), barWidth, 'FaceColor', cols(i,:), 'EdgeColor', 'none');
    
    % draw error bars
    errorbar(xPos, mu(:,i), se(:,i), 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 4);
end
xlabel('Time window')
set(gca,'XTick',1:3)
set(gca,'XTickLabels',timenames)
xlim([.5 3.5])
set(gca,'FontSize',20)
ylabel('RH-LH difference')
title('All stimuli')
ylim(ylims)


% objects/words separately
load('/Users/uqarobi7/Dropbox/PROJECTS/hemiraptor/hemi-central/results/stats_imagecorrs_perobjs.mat')

objs = {'objects' 'words'};
for o=1:2
    subplot(3,2,2+(o-1)*2+1)
    plot(timevect,timevect*0,'LineWidth',0.5,'Color','k','HandleVisibility','off')
    hold on
    diff=struct();
    for m = 1:length(modelnames)
        diff.(modelnames{m}).mu = stats.partialcorr_r.(objs{o}).(modelnames{m}).mu-stats.partialcorr_L.(objs{o}).(modelnames{m}).mu;
        diff.(modelnames{m}).mu_all = stats.partialcorr_r.(objs{o}).(modelnames{m}).mu_all-stats.partialcorr_L.(objs{o}).(modelnames{m}).mu_all;

        plot(timevect,movmean(diff.(modelnames{m}).mu,20),'LineWidth',2,'Color',cols(m,:))
    end
    for t = 1:size(tws,1)
        fill([tws(t,:) fliplr(tws(t,:))],[-.15 -.15 .25 .25],[.7 .7 .7],'FaceAlpha',.2,'EdgeColor','none')
    end
    % legend(modelnames)

    for t = 1:3
        tw = tws(t,:);
        tw = find(timevect>=tw(1)&timevect<=tw(2));
        for m = 1:length(modelnames)
            mdat = mean(diff.(modelnames{m}).mu_all(tw,:),1);
            mu(t,m) = mean(mdat);
            se(t,m) = std(mdat)/sqrt(length(mdat));
        end
        timenames{t} = sprintf('%d to %d ms',tws(t,1),tws(t,2));
    end
    xlabel('Time (ms)')
    set(gca,'FontSize',20)
    ylabel('RH-LH difference')
    title(objs{o})
    ylim(ylims)

    subplot(3,2,2+(o-1)*2+2);
    hold on;
    nGroups = size(tws,1);
    nBars = 6;
    groupWidth = 0.7;
    barWidth = groupWidth / nBars;
    for i = 1:nBars
        % calculate x positions for this bar across both groups
        xPos = (1:nGroups) - groupWidth/2 + (i-0.5) * barWidth;

        % draw bars
        bar(xPos, mu(:,i), barWidth, 'FaceColor', cols(i,:), 'EdgeColor', 'none');

        % draw error bars
        errorbar(xPos, mu(:,i), se(:,i), 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 4);
    end
    xlabel('Time window')
    set(gca,'XTick',1:3)
    set(gca,'XTickLabels',timenames)
    xlim([.5 3.5])
    set(gca,'FontSize',20)
    ylabel('RH-LH difference')
    title(objs{o})
    ylim(ylims)

end