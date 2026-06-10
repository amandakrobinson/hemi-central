
function plot_decoding_central()

%% stuff for plotting
addpath('~/Dropbox/MATLAB/Colormaps/Colormaps (5)/Colormaps')
addpath('~/Dropbox/MATLAB/altmany-export_fig')
colormap viridis
cmap = viridis(3);
cmap(3,:) = [0 0 0];

%% load results
load('results/stats_decoding_central.mat','stats');
timevect = stats.timevect;

%% set variables
clustnames = {'Left hemisphere' 'Right hemisphere' 'RH-LH difference'};
clust = {'L_elecs' 'R_elecs' 'RHLHdiff'};

% plot details
xlims = [-100 800];
yplotlim = [0.48 0.62;-.02 .1];

%% plot central overlaid

% initialise figure
fig=figure(4);clf;
fig.Resize = 'off';
fig.Position=[1 1 800 900];

subplot(2,1,1);

for cl = 1:2

    hold on
    plot(timevect,timevect*0+0.5,'k','HandleVisibility','off')

    col = cmap(cl,:);

    % plot decoding
    dat = stats.(clust{cl});
    mu = dat.mu;
    se = dat.se;

    fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
    plot(timevect,mu,'Color',col,'LineWidth',2)

    % plot onset and peak CIs
    yval = 0.495-cl*.005;

    % plot([dat.onset dat.peak],[yval yval],'.','Color',col,'MarkerSize',12,'HandleVisibility','off')
    % line(dat.onsetci+[-.5 .5],[yval yval],'Color',col,'LineWidth',2,'HandleVisibility','off')
    % line(dat.peakci+[-.5 .5],[yval yval],'Color',col,'LineWidth',2,'HandleVisibility','off')
    % 
    % if cl == 1
    %     text(median(dat.onsetci),yval+.005,'onset  ','HorizontalAlignment','center')
    %     text(median(dat.peakci),yval+.005,'  peak','HorizontalAlignment','center')
    % end
    drawnow

end

set(gca,'FontSize',20)
set(gca,'linewidth',2)
set(gca,'YTick',.5:.02:.65);
set(gca,'TickDir','out');
set(gca,'XTick',0:200:800);

legend({'Left hemisphere' 'Right hemisphere'},'Box','off')
ylim([0.48 0.62]);
xlim([-100 800]);
xlabel('Time (ms)')
ylabel('Decoding accuracy')

bflabels = {'Left hemisphere' 'Right hemisphere'};

% plot bayes factors
for v = 1:2 % new plot for each hemi

    col = cmap(v,:);
    bf = stats.(clust{v}).bf;

    sp=subplot(6,1,3+v); % LH is upper plot
    set(sp,'FontSize',20);
    hold on

    plot(timevect,1+0*timevect,'k-');
    cob = [.5 .5 .5;1 1 1;col];
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
    ylim(10.^(7.5*[-1 1]))
    sp.YTick = 10.^([-7 0 7]);
    xlim([min(timevect) max(timevect)+1])
    ylabel('BF')
    xlabel('Time (ms)')

    % add labels
    bflab = bflabels{v};
    text(500,1e5,bflab,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',col);

end

%% now do difference
% plot bayes factors
sp=subplot(6,1,6);
set(sp,'FontSize',20);
hold on
bf = stats.RHLHdiff.bf;
col = [0 0 0];
plot(timevect,1+0*timevect,'k-');
cob = [.5 .5 .5;1 1 1;col];
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
ylim(10.^(6.5*[-1 1]))
sp.YTick = 10.^([-5 0 5]);
xlim([min(timevect) max(timevect)+1])
ylabel('BF')
xlabel('Time (ms)')

% add labels
text(300,1e4,'Right-left difference','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',col);

% annotate and save
annotation('textbox',[0.06 .89 .1 .1],...
    'String','A','FontSize',30,'LineStyle','none')

annotation('textbox',[0.5 .89 .1 .1],...
    'String','B','FontSize',30,'LineStyle','none')

export_fig('figures/decoding_central_overlaid.png','-transparent','-r300')
