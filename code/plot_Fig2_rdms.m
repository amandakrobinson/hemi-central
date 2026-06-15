%% plot Figure 2: Model RDMs

load('results/RSAmodels.mat','m_cielab','m_sflow','m_sfhigh','m_rect','m_lum', ...
    'm_cat','m_imw','m_edge');

figure(1);clf;
set(gcf,'Position',[1 1 1050 1050],'Resize','off')
cm = viridis;
set(gcf,'Colormap',cm)

cmap = tab10(8);
medcat = [2.5:4:24 25.5:2:36];
stimorder = {'fish','bird','face','boats','tree','tools'};

fs = 14;

subplot(3,3,1);
imagesc(m_edge,[0.9 1]);
axis square
set(gca,'FontSize',fs)
set(gca,'YDir','normal');
set(gca,'XTick',[])
set(gca,'YTick',medcat)
set(gca,'YTickLabels',[stimorder stimorder])
title('Edge','Color',cmap(1,:),'FontSize',20)

subplot(3,3,2)
imagesc(m_cielab);
set(gca,'YDir','normal');
set(gca,'FontSize',fs)
title('Colour','Color',cmap(2,:),'FontSize',20)
axis square
axis off

subplot(3,3,3)
imagesc(m_lum)
set(gca,'YDir','normal');
set(gca,'FontSize',fs)
title('Mean luminance','Color',cmap(3,:),'FontSize',20)
axis square
axis off

subplot(3,3,4)
imagesc(m_sflow)
set(gca,'FontSize',fs)
title('Low spatial frequency','Color',cmap(4,:),'FontSize',20)
axis square
set(gca,'YDir','normal');
set(gca,'XTick',[])
set(gca,'YTick',medcat)
set(gca,'YTickLabels',[stimorder stimorder])

subplot(3,3,5)
imagesc(m_sfhigh)
set(gca,'YDir','normal');
axis square
set(gca,'FontSize',fs)
axis off
title('High spatial frequency','Color',cmap(5,:),'FontSize',20)

subplot(3,3,6)
imagesc(m_rect,[0 3])
set(gca,'YDir','normal');
axis square
set(gca,'FontSize',fs)
title('Rectilinearity','Color',cmap(6,:),'FontSize',20)
axis off

subplot(3,3,7)
imagesc(m_cat);
set(gca,'FontSize',fs)
title('Concept category','Color',cmap(7,:),'FontSize',20)
axis square
set(gca,'YDir','normal');
set(gca,'XTick',medcat)
set(gca,'XTickLabels',[stimorder stimorder])
set(gca,'YTick',medcat)
set(gca,'YTickLabels',[stimorder stimorder])
set(gca,'XTickLabelRotation',90)

subplot(3,3,8)
imagesc(m_imw)
set(gca,'YDir','normal');
set(gca,'FontSize',fs)
title('Image/word','Color',cmap(8,:),'FontSize',20)
axis square
set(gca,'XTick',medcat)
set(gca,'XTickLabels',[stimorder stimorder])
set(gca,'XTickLabelRotation',90)
set(gca,'YTick',[])

% add colourbar without resizing plot
ax = gca;
orig_pos = ax.Position;
cb = colorbar('southoutside');
% cb.Position = cb.Position + 1e-10;
cb.Position(2) = .02;
cb.Ticks = cb.Limits;               % Set specific tick locations
cb.TickLabels = {'Similar', 'Dissimilar'}; % Change text for those ticks
ax.Position = orig_pos;


% correlations between
mask = ones(size(m_cielab));
mask = tril(mask,-1);

mods = [m_edge(mask==1) m_cielab(mask==1) m_lum(mask==1) m_sflow(mask==1) m_sfhigh(mask==1) m_rect(mask==1) ...
    m_cat(mask==1) m_imw(mask==1)];
labels = {'Edge' 'Colour',...
    'Mean luminance','Low spatial frequency', ...
    'High spatial frequency','Rectilinearity',...
    'Concept','Image/Word'};

dists = 1-corr(mods,'type','spearman');

a=mdscale(dists,2,'Criterion','metricstress');

% labels for left plots
for np = 1:3
    yadjust = (np-1)*.3;
    annotation(gcf,'textarrow',...
        [0.06 0.06], [0.785-yadjust 0.8-yadjust],...
        'String','Images', 'HeadStyle', 'none', 'LineStyle', 'none',...
        'FontSize',20, 'color','k',...
        'FontWeight','bold', 'TextRotation',90);
    annotation(gcf,'textarrow',...
        [0.06 0.06], [0.888-yadjust 0.95-yadjust],...
        'String','Words', 'HeadStyle', 'none', 'LineStyle', 'none',...
        'FontSize',20, 'color','k',...
        'FontWeight','bold', 'TextRotation',90);
end

% labels for bottom plots
for np = 1:2
    xadjust = (np-1)*.28;
    annotation(gcf,'textarrow',...
        [0.24+xadjust 0.3+xadjust], [0.06 0.06],...
        'String','Images', 'HeadStyle', 'none', 'LineStyle', 'none',...
        'FontSize',20, 'color','k',...
        'FontWeight','bold');
    annotation(gcf,'textarrow',...
        [0.34+xadjust 0.41+xadjust], [0.06 0.06],...
        'String','Words', 'HeadStyle', 'none', 'LineStyle', 'none',...
        'FontSize',20, 'color','k',...
        'FontWeight','bold');
end

% plot model mds
sp=subplot(3,3,9); hold on;
set(gca,'FontSize',18)

xs = [.05 0 0 0 0 -.05 0 0];
ys = [0 -.08 -.08 -.12 -.12 .08 -.08 -.08];
align = {'left' 'center' 'center' 'center' 'center' 'left' 'center' 'center'};
labels{4} = {'Low spatial','frequency'};
labels{5} = {'High spatial','frequency'};
for l = 1:length(labels)
    scatter(a(l,1),a(l,2),50,cmap(l,:),'filled')
    text(a(l,1)+xs(l),a(l,2)+ys(l),labels{l},'HorizontalAlignment',align{l},'FontSize',12)
end
xlim([-.7 1]);
ylim([-.6 .8]);
sp.XTick = [];
sp.YTick = [];
xlabel('Dimension 1')
ylabel('Dimension 2')
title('Model relationships','FontSize',20)

%% save
fn = './figures/Figure2_models';

% export PNG
set(gcf, 'Renderer', 'painters')
print(gcf, '-dpng', '-r500', fn)
im = imread([fn '.png']);
[i,j] = find(mean(im,3) < 255); margin = 2;
imwrite(im(min(i)-margin:max(i)+margin, min(j)-margin:max(j)+margin, :), [fn '.png'], 'png')

%% correlations between models

% in order of 'labels'
[r,p] = corr(mods,'type','spearman');
models = {'Edge' 'Colour',...
    'Mean luminance','Low spatial frequency', ...
    'High spatial frequency','Rectilinearity',...
    'Concept','Image/Word'};


tablecorrs = array2table(round(r,2),'VariableNames',validlabels','RowNames',models');
tablepvals = array2table(p,'VariableNames',validlabels','RowNames',models');
