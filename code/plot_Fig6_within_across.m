
%% load data for correlating 
load('results/stats_withinacross.mat','stats','timevect');

%% plot within and across timegens

figure(2);clf
set(gcf,'Position',[674 92 1200 1000],'Resize','off')

cmaplabel = viridis(3);
cmaplabel(3,:) = [.5 .5 .5];

cmap = brewermap(50,'Oranges');

cmap(1,:) = [1 1 1];

timv = timevect;
timerv = fliplr(timv);
xlab = 0:100:500;
ylab = 500:-100:0;

ts = [0 .3]; % scales

conds = {'L_elecs' 'R_elecs' 'across'};
labels = {'Within LH' 'Within RH' 'Across LH:RH'};

for c = 1:length(conds)

    a=subplot(3,3,c);hold on

    dat = stats.(conds{c});
    mu = dat.mu;
    bf = dat.bf;
    mu(bf<10) = 0;

    %%% plot time gen decoding
    ax=imagesc(timv,timv,(mu),ts);
    hold on
    plot(timv,timv,'Color',[.5 .5 .5])
    colormap(gca,cmap)

    set(gca,'YDir','normal')
    cb=colorbar();
    y=ylabel(cb,'Correlation','FontSize',14,'Rotation',270);
    y.Position(1) = cb.Position(1)+2;

    cb.Ticks = [ts(1) ts(2)];
    xlim([0 500])
    ylim([0 500])
    set(gca,'XTick',[0:100:500 max(timv)],'XTickLabel',num2cell(xlab));
    set(gca,'YTick',[0:100:500 max(timv)],'YTickLabel',num2cell(xlab));

    if c < 3
        xlabel('Time (ms) - odd sequences');
        ylabel('Time (ms) - even sequences');
    else
        xlabel('Left hemisphere time (ms)');
        ylabel('Right hemisphere time (ms)');
    end

    axis square
    title(labels{c},'Color',cmaplabel(c,:))
    set(gca,'Fontsize',14)

end

%% plot within-across differences

diffts = [-.2 .2];
cmap = brewermap(50,'PuOr');
cmap = flipud(cmap);

difflabels = {'Within LH-across difference'; 'Within RH-across difference'};

for c = 1:2

    a=subplot(3,3,3+c);hold on

    dat = stats.withinVacross.(conds{c});
    mu = dat.mu;%reshape(dat.mu,length(timv),length(timv));
    bf = dat.bf;%reshape(dat.bf,length(timv),length(timv));
    mu(bf<10) = 0;

    %%% plot time gen decoding
    ax=imagesc(timv,timv,(mu),diffts);
    hold on
    plot(timv,timv,'Color',[.5 .5 .5])
    colormap(gca,cmap)

    set(gca,'YDir','normal')
    cb=colorbar();
    y=ylabel(cb,'Correlation difference','FontSize',14,'Rotation',270);
    y.Position(1) = cb.Position(1)+2;

    cb.Ticks = [diffts(1) diffts(2)];
    xlim([0 500])
    ylim([0 500])
    set(gca,'XTick',[0:100:500 max(timv)],'XTickLabel',num2cell(xlab));
    set(gca,'YTick',[0:100:500 max(timv)],'YTickLabel',num2cell(xlab));

    xlabel('Time (ms)');
    ylabel('Time (ms)');

    axis square
    title(difflabels{c})
    set(gca,'Fontsize',14)

end


%% plot RH-LH differences

diffts = [-.2 .2];
cmap = brewermap(50,'PuOr');
cmap = flipud(cmap);

difflabels = {'RH-LH difference'};

a=subplot(3,3,6);hold on

dat = stats.LHRHdiff;
mu = dat.mu;
bf = dat.bf;
mu(bf<10) = 0;

%%% plot RH-LH diff
ax=imagesc(timv,timv,(mu),diffts);
hold on
plot(timv,timv,'Color',[.5 .5 .5])
colormap(gca,cmap)

set(gca,'YDir','normal')
cb=colorbar();
y=ylabel(cb,'Correlation difference','FontSize',14,'Rotation',270);
y.Position(1) = cb.Position(1)+2;

cb.Ticks = [diffts(1) diffts(2)];
xlim([0 500])
ylim([0 500])
set(gca,'XTick',[0:100:500 max(timv)],'XTickLabel',num2cell(xlab));
set(gca,'YTick',[0:100:500 max(timv)],'YTickLabel',num2cell(xlab));

xlabel('Time (ms)');
ylabel('Time (ms)');

axis square
title(difflabels{1})
set(gca,'Fontsize',14)

%% plot diagonal 
% line plot for within/across for withinLH, within RH and across

% plot details
xlims = [-100 800];
yplotlim = [-.04 0.4];

comps = {'Within LH' 'Within RH' 'Across LH:RH'};
bfcomps = {'Within versus across: LH' 'Within versus across: RH'};

% plot correlations
sub=subplot(6,2,9);cla
sub.Position(2) = sub.Position(2)-.019;
set(gca,'FontSize',14)
set(gca,'linewidth',2)
sub.YTick = 0:.1:1;
sub.TickDir = 'out';
sub.XTick = 0:200:1000;

hold on
plot(timevect,timevect*0,'Color',[.7 .7 .7],'HandleVisibility','off')
plot(yplotlim*0,yplotlim,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

% get dat
cont = {'L_elecs' 'R_elecs' 'across'};
cmap = cmaplabel;

for x = 1:3
    dat = stats.(cont{x});

    col = cmap(x,:);

    mu = diag(dat.mu)';
    se = diag(dat.se)';

    fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
    plot(timevect,mu,'Color',col,'LineWidth',2)
    drawnow
end

ylim(yplotlim);
xlim(xlims);
% xlabel('Time (ms)')
ylabel('Correlation')

l=legend(comps,'Box','off');
% l.Position(1) = l.Position(1)+.005;


% plot different BF
for c = 1:2
    sp=subplot(12,2,(c-1)*2+21);
    sp.Position(2) = sp.Position(2)-.019;

    % sp.Position(3) = sp.Position(3)+.085;
    % sp.Position(4) = .05;

    set(sp,'FontSize',14);
    hold on

    dat = stats.withinVacross.(cont{c});

    mu = diag(dat.mu);
    bf = diag(dat.bf);

    idxmu = [mu<0,mu>0]';

    plot(timevect,1+0*timevect,'k-');
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    for m=1:2
        col = cmap(c,:);
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
    ylim([10.^-4 10.^7])
    sp.YTick = 10.^([-4 0 4]);
    sp.XTick = 0:200:800;
    sp.TickDir = 'out';
    xlim([min(timevect) max(timevect)+1])
    if c == 2
        xlabel('Time (ms)')
        yl=ylabel('BF_d_i_f_f');
        set(yl,'Units','normalized')
        yl.Position(2) = yl.Position(2)+.9;
    else
        sp.XTickLabel = [];
    end
end

%% plot correlations by time offset

offtimes = -150:150;

trange = find(timevect>=0&timevect<=500);

offdat=struct();
for p = 1:20
    x=0;
    for ot = offtimes
        x=x+1;
        offdat.muL(x,p) = mean(diag(stats.L_elecs.mu_all(trange,trange,p),ot));
        offdat.muR(x,p) = mean(diag(stats.R_elecs.mu_all(trange,trange,p),ot));
        offdat.muA(x,p) = mean(diag(stats.across.mu_all(trange,trange,p),ot));
    end
end


% stats
cont = {'muL' 'muR' 'muA'};

for x = 1:2
    dat = offdat.(cont{x})-offdat.muA; %within minus across

    mu = mean(dat,2)';
    offdat.bf.(cont{x}) = bayesfactor_R_wrapper(dat,...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
end


%% plot
cmap = cmaplabel;
cont = {'muL' 'muR' 'muA'};

% plot correlations
sub=subplot(6,2,10);cla
sub.Position(2) = sub.Position(2)-.019;
% sub.Position(3) = sub.Position(3)+.085;

set(gca,'FontSize',14)
set(gca,'linewidth',2)
sub.YTick = 0:.1:2;
sub.TickDir = 'out';
sub.XTick = -150:50:150;

yplotlim = [-.005 .2];

hold on
plot(offtimes,offtimes*0,'Color',[.7 .7 .7],'HandleVisibility','off')
plot(yplotlim*0,yplotlim,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms


for x = 1:3
    dat = offdat.(cont{x});

    col = cmap(x,:);

    mu = mean(dat,2)';
    se = std(dat,[],2)'/sqrt(size(dat,2));

    fill([offtimes fliplr(offtimes)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
    plot(offtimes,mu,'Color',col,'LineWidth',2)
    drawnow

end

ylim(yplotlim);
xlim([min(offtimes) max(offtimes)]);
% xlabel('Time delay (ms)')
ylabel('Mean correlation')

l=legend(comps,'Box','off');
% l.Position(1) = l.Position(1)+.005;

%% plot different BF
cont = {'muL' 'muR' 'muA'};

for c = 1:2
    sp=subplot(12,2,(c-1)*2+22);

    % if c == 1
    %     sp.Position(2) = sp.Position(2)-.018;
    % else 
    %     sp.Position(2) = sp.Position(2)+.01;
    % end
    % 
    sp.Position(2) = sp.Position(2)-.019;
    % sp.Position(4) = .05;

    set(sp,'FontSize',14);
    hold on

    bf = offdat.bf.(cont{c});

    plot(offtimes,1+0*offtimes,'k-');
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    col = cmap(c,:);
    cob = [.5 .5 .5;1 1 1;col];
    for i=1:3
        x = offtimes(idx(i,:));
        y = bf(idx(i,:));
        if ~isempty(x)
            stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
        end
    end

    sp.YScale='log';
    ylim([10.^-4 10.^7])
    sp.YTick = 10.^([-4 0 4]);
    sp.XTick = -150:50:150;
    sp.TickDir = 'out';

    xlim([min(offtimes) max(offtimes)])
    if c == 2
        xlabel('Time delay (ms)')
        yl=ylabel('BF_d_i_f_f');
        set(yl,'Units','normalized')
        yl.Position(2) = yl.Position(2)+0.9;
    else
        sp.XTickLabel = {};
    end
end

%% annotate
annotation('textbox',[0.07 .85 .1 .1],...
    'String','A','FontSize',30,'LineStyle','none')

annotation('textbox',[0.07 .56 .1 .1],...
    'String','B','FontSize',30,'LineStyle','none')

annotation('textbox',[0.07 .27 .1 .1],...
    'String','C','FontSize',30,'LineStyle','none')

annotation('textbox',[0.5 .27 .1 .1],...
    'String','D','FontSize',30,'LineStyle','none')

%% save
fn = './figures/Figure6_withinacross';
set(gcf, 'Renderer', 'painters')
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

