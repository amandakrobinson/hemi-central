% plot RDMs and MDS for left/right hemi neural representations at 120-150ms

%% get neural rdms

load('results/rdms.mat','alldat')
timevect = alldat.time;
clust = {'L_elecs' 'R_elecs'};
clustnames = {'Left hemisphere' 'Right hemisphere'};
load('results/imageorder.mat','imnames')

medcat = [2.5:4:24 25.5:2:36];
stimorder = {'fish','bird','face','boats','tree','tools'};

%% plot MDS

timetoplot = find(timevect>119&timevect<151);

rng(7302)

f=figure(1);clf;
f.Position = [f.Position(1:2) 950 800];

meanrdm = (alldat.(clust{1})+alldat.(clust{2}))/2;
meanrdm = mean(mean(meanrdm(:,timetoplot,:),2),3);
Xstart = mdscale(squareform(meanrdm),2,'Start','random','Replicates',50);

for h = 1:length(clust)

    % RDM
    a = subplot(2,2,h);
    a.FontSize = 18;
    rdm = mean(mean(alldat.(clust{h})(:,timetoplot,:),2),3);
    imagesc(squareform(rdm),[0.5 0.65]);axis square
    a.YDir='normal';
    colorbar
    colormap viridis
    set(a,'Fontsize',18)
    set(gca,'YTick',medcat)
    set(gca,'YTickLabels',[stimorder stimorder])
    set(gca,'XTick',medcat)
    set(gca,'XTickLabels',[stimorder stimorder])
    set(gca,'XTickLabelRotation',90)
    t=title(sprintf('%s: %d-%d ms',clustnames{h},min(timevect(timetoplot)),max(timevect(timetoplot))),...
        'FontSize',22);
    t.Position(2) = t.Position(2)+1;

    % embedding
    fprintf('create embedding\n')
    Y = squareform(rdm);

    if h == 1
        X1 = mdscale(Y,2,'Start',Xstart);
    else
        X1 = mdscale(Y,2,'Start',X1);
        X1(:,1) = -1*X1(:,1); % reverse X axis so looks visually similar
    end

    a = subplot(2,2,h+2);

    a.FontSize = 18;
    a.LineWidth = 2;
    hold on
    X=X1;
    aw=range(X(:))*.06;
    imrange = 1:128;
    plot(X(:,1),X(:,2),'.')
    for i=1:length(imnames)
        imname = dir(sprintf('stimuli_rendered/stim%03i*',i));
        IM = imread(sprintf('stimuli_rendered/%s',imname.name));
        IM = flipud(IM);
        if size(IM,3)==1
            IM = repmat(IM,1,1,3);
        end
        
        image('XData',X(i,1)+aw*[-1 1],'YData',X(i,2)+aw*[-1 1],'CData',IM(imrange,imrange,1:3))
    end
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    axis equal
    % axis square

    a.XTick = [];
    a.YTick = [];
    ylim([-.38 .41])
    xlim([-.41 .38])
end
drawnow

% after the loop ends, adjust positions of bottom subplots
for h = 1:length(clust)
    ax = subplot(2,2,h+2);
    ax.Position(1) = ax.Position(1) - 0.02;
end


% labels for left plots
annotation(gcf,'textarrow',...
    [0.06 0.06], [0.74 0.74],...
    'String','Images', 'HeadStyle', 'none', 'LineStyle', 'none',...
    'FontSize',20, 'color','k',...
    'FontWeight','bold', 'TextRotation',90);
annotation(gcf,'textarrow',...
    [0.06 0.06], [0.9 .9],...
    'String','Words', 'HeadStyle', 'none', 'LineStyle', 'none',...
    'FontSize',20, 'color','k',...
    'FontWeight','bold', 'TextRotation',90);


% labels for bottom plots
for np = 1:2
    xadjust = (np-1)*.44;
    annotation(gcf,'textarrow',...
        [0.27+xadjust 0.27+xadjust], [0.5 0.5],...
        'String','Images', 'HeadStyle', 'none', 'LineStyle', 'none',...
        'FontSize',20, 'color','k',...
        'FontWeight','bold');
    annotation(gcf,'textarrow',...
        [0.38+xadjust 0.38+xadjust], [0.5 0.5],...
        'String','Words', 'HeadStyle', 'none', 'LineStyle', 'none',...
        'FontSize',20, 'color','k',...
        'FontWeight','bold');
end


%% save
fn = 'figures/Figure4_neuralMDS';

% export PNG
set(gcf, 'Renderer', 'painters')
print(gcf, '-dpng', '-r500', fn)
im = imread([fn '.png']);
[i,j] = find(mean(im,3) < 255); margin = 2;
imwrite(im(min(i)-margin:max(i)+margin, min(j)-margin:max(j)+margin, :), [fn '.png'], 'png')