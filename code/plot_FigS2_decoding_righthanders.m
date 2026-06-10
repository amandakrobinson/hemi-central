
function plot_FigS2_decoding_righthanders()

%% stuff for plotting
cmap = viridis(256);
co3 = cmap(round(linspace(1, 256, 3)), :);
co3(3,:) = [0 0 0];
cmap = co3;

%% load results
load('results/stats_decoding_central.mat','stats');
timevect = stats.timevect;

%% get participant details
participants = readtable('../hemifield-eeg/data/participants.tsv',"FileType","text",'Delimiter', '\t');
righthanders = participants.Total>50;

%% set variables
clust = {'L_elecs' 'R_elecs' 'RHLHdiff'};

%% plot central overlaid

% initialise figure
fig=figure(4);clf;
fig.Resize = 'off';
fig.Position=[1 1 800 800];

subplot(2,1,1);

for cl = 1:2

    hold on
    plot(timevect,timevect*0+0.5,'k','HandleVisibility','off')

    col = cmap(cl,:);

    % plot decoding
    dat = stats.(clust{cl});
    x = dat.mu_all(righthanders,:);
    mu = mean(x,1);
    se = std(x)./sqrt(size(x,1));

    fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
    plot(timevect,mu,'Color',col,'LineWidth',2)

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

%%
bflabels = {'Left hemisphere' 'Right hemisphere'};

% plot bayes factors
for v = 1:2 % new plot for each hemi

    col = cmap(v,:);
    dat = stats.(clust{v});
    % calculate bayesfactors
    bf = bayesfactor_R_wrapper(dat.mu_all(righthanders,:)',...
        'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

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
dat = stats.RHLHdiff;
% calculate bayesfactors
bf = bayesfactor_R_wrapper(dat.mu_all(righthanders,:)',...
    'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

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

%% save
fn = './figures/FigureS2_decoding_righthanders';
set(gcf, 'Renderer', 'painters')
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
