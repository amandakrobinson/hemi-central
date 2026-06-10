% plot Fig5 image RSA

%% load results
load('results/stats_imagecorrs.mat','stats','modelnames');
timevect = stats.timevect;
features = fieldnames(stats.partialcorr_L);

% set variables
clustnames = {'Left hemisphere' 'Right hemisphere' 'RH-LH difference'};

%% plotting each model on a separate plot
% stuff for plotting labels
cmap = tab10(6);
cmap(7,:) = [0 0 0];

% always sample from the full viridis range consistently
vir = viridis(256);
co3 = vir(round(linspace(1, 256, 3)), :);

% set variables
clust = {'partialcorr_L' 'partialcorr_r'};
labels = {'Edge' 'Colour' 'Luminance' 'Low spatial frequency' 'High spatial frequency' 'Rectilinearity'};

% plot details
xlims = [-100 800];
yplotlim = [-.1 0.35];

spl = [1 2 3 10 11 12];

% initialise figure
fig=figure(5);clf;
fig.Resize = 'off';
fig.Position=[1 1 1400 800];

for a = 1:length(modelnames)

    subplot(6,3,[0 3]+spl(a))
    hold on
    plot(timevect,timevect*0,'k','HandleVisibility','off')

    for c = 1:length(clust)

        col = co3(c,:);

        % plot correlation
        dat = stats.(clust{c}).(modelnames{a});
        mu = dat.mu';
        se = dat.se';

        fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
        plot(timevect,mu,'Color',col,'LineWidth',2)

    end
    legend({'LH' 'RH'},'Box','off')
    title(labels{a},'Color',cmap(a,:))
    set(gca,'YTick',[0:.1:.4])
    set(gca,'XTickLabel',{})
    xlim(xlims)
    ylim(yplotlim)
    set(gca,'FontSize',18)
    ylabel('Partial correlation')
    % xlabel('Time (ms)')

    %% plot different BF
    sp=subplot(6,3,6+spl(a));
    set(sp,'FontSize',18);
    hold on

    col = [0 0 0];

    bf = stats.diff_partialcorr.(modelnames{a}).bf;
    mu = stats.diff_partialcorr.(modelnames{a}).mu;
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
    ylim(10.^(5*[-1 1]))
    sp.YTick = 10.^([-5 0 5]);
    sp.XTick = 0:200:800;
    % sp.XTickLabel = [];
    xlim([min(timevect) max(timevect)+1])
    ylabel('BF_d_i_f_f')
    xlabel('Time (ms)')

end

%%
for a= 1:length(modelnames)
    sp=subplot(6,3,6+spl(a));
    sp.Position(2) = sp.Position(2)+.03;
    % sp.Position(4) = sp.Position(4)-.04;
end

%% save
fn = './figures/Figure5_imageRSA';
set(gcf, 'Renderer', 'painters')
print(gcf, '-dpng', '-r500', fn)
im = imread([fn '.png']);
[i,j] = find(mean(im,3) < 255); margin = 2;
imwrite(im(min(i)-margin:max(i)+margin, min(j)-margin:max(j)+margin, :), [fn '.png'], 'png')
