
% plot Figures : mean decoding accuracy per object class for central stimuli

function plot_decoding_perobject()

%% stuff for plotting
addpath('~/Dropbox/MATLAB/Colormaps/Colormaps (5)/Colormaps')
addpath('~/Dropbox/MATLAB/altmany-export_fig')
co3 = viridis(3);
co3(3,:) = [0 0 0];

%% load results
load('results/stats_decoding_perobs.mat','stats');
timevect = stats.timevect;

%% set variables
combos = {'Face' 'Word' 'Object' 'FaceObject' 'FaceWord' 'WordObject'};
clustnames = {'Left hemisphere' 'Right hemisphere' 'Right-Left difference'};

% plot details
xlims = [-100 800];
yplotlim = [0.48 0.665];
chance = [0.5 0.5 0];

%% plot paired decoding results

% initialise figure
fig=figure(1);clf;
fig.Resize = 'off';
fig.Position=[1 1 1800 1000];

plotnr = [14 27 1 13 26 25];
plotnr = [2 3 1 13 15 14];


for o = 1:length(combos)

    st = plotnr(o);
    sub=subplot(12,3,[st st+3]);

    %% plot per condition
    dat_p = stats.(combos{o});
    clust = fieldnames(dat_p); % LH, RH, LH-RH diff

    for cl = 1:2

        % get data to plot
        dat_cond = dat_p.(clust{cl});

        set(gca,'FontSize',18)
        set(gca,'linewidth',2)
        sub.YTick = .5:.05:.65;

        sub.TickDir = 'out';
        sub.XTick = 0:200:1000;
        sub.XTickLabel = {};

        hold on
        plot(timevect,timevect*0+chance(cl),'k','HandleVisibility','off')
        plot(yplotlim*0,yplotlim,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

        col = co3(cl,:);

        mu = dat_cond.mu;
        se = dat_cond.se;

        fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
        plot(timevect,mu,'Color',col,'LineWidth',2);

    end

    ylim(yplotlim);
    xlim(xlims);
    title(combos{o},'FontSize',25)
    % xlabel('Time (ms)')
    ylabel('Accuracy','FontSize',18)

    if st==1
        l=legend(clustnames(1:2),'Box','off','FontSize',18);
    end
end
% set(l,'Position',[.4 .82 .1 .1])

%% plot bayes factors
% initialise figure
% fig=figure(2);clf;
% fig.Resize = 'off';
% fig.Position=[1 1 2000 900];

for o = 1:length(combos)

    sp=subplot(12,3,plotnr(o)+6);
    sp.Position(2) = sp.Position(2)+.02;
    sp.Position(4) = sp.Position(4)-.02;

    set(gca,'FontSize',18)

    hold on
    %% plot per condition
    dat_p = stats.(combos{o});

    % get data to plot
    bf = dat_p.LHRHdiff.bf;

    plot(timevect,1+0*timevect,'k-');
    cob = [.5 .5 .5;1 1 1;0 0 0];
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    for i=1:3
        x = timevect(idx(i,:));
        y = bf(idx(i,:));
        if ~isempty(x)
            stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
        end
    end
    sp.YScale='log';
    ylim(10.^(4*[-1 1]))
    sp.YTick = 10.^([-4 0 4]);
    sp.YAxis.FontSize = 12
    xlim([min(timevect) max(timevect)+1])
    ylabel('BF','FontSize',18)
    xlabel('Time (ms)')
end

%% save
fn = './figures/Figure_decoding_perobject';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');



end


