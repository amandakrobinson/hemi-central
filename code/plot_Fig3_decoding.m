
% plot Figures : mean decoding accuracy per object class for central stimuli

function plot_Fig3_decoding()

%% stuff for plotting

% always sample from the full viridis range consistently
cmap = viridis(256);
co3 = cmap(round(linspace(1, 256, 3)), :);
co3(3,:) = [0 0 0];

%% PLOT DECODING FOR MEAN OF ALL 630 stimulus pairs
% load results
load('results/stats_decoding_central.mat','stats');
timevect = stats.timevect;

% set variables
clust = {'L_elecs' 'R_elecs' 'RHLHdiff'};
bflabels = {'Left hemisphere' 'Right hemisphere'};

% plot central overlaid

% initialise figure
fig=figure(1);clf;
fig.Units = 'pixels';
fig.Position=[1 1 1500 900];
fig.Resize = 'off';

% plot details
xlims = [-100 800];
yplotlim = [0.48 0.62];

sp1=subplot(2,2,1);
sp1.Position(1) = sp1.Position(1)-.05;
sp1.Position(3) = sp1.Position(3)+.05;   
for cl = 1:2

    hold on
    plot(timevect,timevect*0+0.5,'k','HandleVisibility','off')

    col = co3(cl,:);

    % plot decoding
    dat = stats.(clust{cl});
    mu = dat.mu;
    se = dat.se;

    fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
    plot(timevect,mu,'Color',col,'LineWidth',2)

    drawnow

end

set(gca,'FontSize',18)
set(gca,'linewidth',2)
set(gca,'YTick',.5:.02:.65);
set(gca,'TickDir','out');
set(gca,'XTick',0:100:800);

text(98,.523,bflabels{1},'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',co3(1,:));
text(280,.57,bflabels{2},'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',co3(2,:));

ylim(yplotlim);
xlim(xlims);
xlabel('Time (ms)')
ylabel('Decoding accuracy')


% plot bayes factors
for v = 1:2 % new plot for each hemi

    col = co3(v,:);
    bf = stats.(clust{v}).bf;

    sp=subplot(6,2,(v-1)*2+7); % LH is upper plot
    sp.Position(1) = sp.Position(1)-.05;
    sp.Position(3) = sp.Position(3)+.05;    
    sp.Position(2) = sp.Position(2)+v*.02;
    sp.Position(4) = sp.Position(4)-.015;
    set(sp,'FontSize',18);
    hold on

    plot(timevect,1+0*timevect,'k-');
    cob = [.5 .5 .5;1 1 1;col];
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    for i=1:3
        x = timevect(idx(i,:));
        y = bf(idx(i,:));
        if ~isempty(x)
            s=stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',2,'MarkerFaceColor',cob(i,:),'LineWidth',.02,'Clipping','on');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
        end
    end
    sp.YScale='log';
    ylim(10.^(7.5*[-1 1]))
    sp.YTick = 10.^([-5 0 5]);
    sp.XTick = 0:100:800;
    sp.XTickLabel = [];
    xlim([min(timevect) max(timevect)+1])

    % add labels
    bflab = bflabels{v};
    text(100,5e-4,bflab,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',col);
    sp.YAxis.FontSize = 16;
    if v == 2
        ylabel('BF','FontSize',22)
    end
end

% now do difference bayes factors
sp=subplot(6,2,11);
sp.Position(1) = sp.Position(1)-.05;
sp.Position(3) = sp.Position(3)+.05;
sp.Position(2) = sp.Position(2)+3*.02;
sp.Position(4) = sp.Position(4)-.015;
hold on
bf = stats.RHLHdiff.bf;
mu = stats.RHLHdiff.mu';
idxmu = [mu<0,mu>0]';
plot(timevect,1+0*timevect,'k-');
idx = [bf<1/10,1/10<bf & bf<10,bf>10]';

for m=1:2
    col = co3(m,:);
    cob = [.5 .5 .5;1 1 1;col];
    for i=1:3
        x = timevect(idx(i,:)&idxmu(m,:));
        y = bf(idx(i,:)&idxmu(m,:));
        if ~isempty(x)
            stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
        end
    end
end
sp.YScale='log';
ylim(10.^(6.5*[-1 1]))
sp.YTick = 10.^([-5 0 5]);    

xlabel('Time (ms)')
sp.XTick = -100:100:800;    
set(sp,'FontSize',18);

sp.YAxis.FontSize = 16;
xlim([min(timevect) max(timevect)+1])

% add labels
text(100,5e-4,'Right-left difference','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',[0 0 0]);

%% PLOT DECODING PER CATEGORY 
% load results
load('results/stats_decoding_perobs.mat','stats');

% set variables
combos = {'Object' 'Word' 'Face' 'WordObject' 'FaceObject' 'FaceWord' };
combonames = {'Objects' 'Words' 'Faces' 'Word vs Object' 'Face vs Object' 'Face vs Word'};

% plot details
xlims = [-100 800];
yplotlim = [0.48 0.665];
chance = [0.5 0.5 0];

% plot paired decoding results
plotnr = [3 4 15 16 27 28];

labelposx = [500 550 550 280 300 340];

for o = 1:length(combos)

    st = plotnr(o);
    sub=subplot(9,4,[st st+4]);
    % if o < 4
    %     sub.Position(2) = sub.Position(2)+.02;
    % end

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
    ylabel('Accuracy','FontSize',18)
% add labels
text(labelposx(o),.65,combonames{o},'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',18,'Color','k','FontWeight','bold');

end

%% bfs
for o = 1:length(combos)

    sp=subplot(9,4,plotnr(o)+8);
    set(gca,'FontSize',18)

    hold on
    %% plot per condition
    dat_p = stats.(combos{o});

    % get data to plot
    bf = dat_p.LHRHdiff.bf; % this is actually RH-LH
    mu = dat_p.LHRHdiff.mu'; % this is actually RH-LH
    idxmu = [mu<0,mu>0]';

    plot(timevect,1+0*timevect,'k-');
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    for m=1:2
        col = co3(m,:);
        cob = [.5 .5 .5;1 1 1;col];
        for i=1:3
            x = timevect(idx(i,:)&idxmu(m,:));
            y = bf(idx(i,:)&idxmu(m,:));
            if ~isempty(x)
                stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
                plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
            end
        end
    end
    sp.YScale='log';
    ylim(10.^(4*[-1 1]))
    sp.YTick = 10.^([-4 0 4]);
    sp.YAxis.FontSize = 12;
    xlim([min(timevect) max(timevect)+1])
    ylabel('BF','FontSize',18)
    xlabel('Time (ms)')
    % sp.Position(4) = .03;
end
%% annotate
annotation('textbox',[0.025 .86 .1 .1],...
    'String','A','FontSize',35,'LineStyle','none')

annotation('textbox',[0.48 .86 .1 .1],...
    'String','B','FontSize',35,'LineStyle','none')

%% save
fn = './figures/Figure3_decoding_all';
set(gcf, 'Renderer', 'painters')
print(gcf, '-dpng', '-r500', fn)
im = imread([fn '.png']);
[i,j] = find(mean(im,3) < 255); margin = 2;
imwrite(im(min(i)-margin:max(i)+margin, min(j)-margin:max(j)+margin, :), [fn '.png'], 'png')

end


