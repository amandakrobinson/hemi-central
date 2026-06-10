
load('results/stats_imagecorrs_perobjs.mat','stats','dat','timevect','mods','modelnames','objconds');

%% stuff for plotting
cmap = viridis(256);
co3 = cmap(round(linspace(1, 256, 3)), :);
co3(3,:) = [0 0 0];

%% PLOT IMAGE STATISTICS CORRELATIONS

% set variables
clust = {'partialcorr_L' 'partialcorr_r'};
labels = {'Edge' 'Colour' 'Luminance' 'Low spatial frequency' 'High spatial frequency' 'Rectilinearity'};

% plot details
xlims = [-100 800];
yplotlim = [-.1 0.35];

spl = [1 2 3 19:21];

for o = 1:length(objconds)

    % initialise figure
    fig=figure(o);clf;
    fig.Resize = 'off';
    fig.Position=[1 1 1400 1000];

    for a = 1:length(modelnames)

        subplot(12,3,[spl(a) spl(a)+3])
        hold on
        plot(timevect,timevect*0,'k','HandleVisibility','off')

        for c = 1:length(clust)

            col = co3(c,:);

            % plot correlation
            dat = stats.(clust{c}).(objconds{o}).(modelnames{a});
            mu = dat.mu';
            se = dat.se';

            fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
            plot(timevect,mu,'Color',col,'LineWidth',2)

        end
        legend({'LH' 'RH'})
        title(labels{a})
        xlim(xlims)
        ylim(yplotlim)
        set(gca,'FontSize',18)
        ylabel('Partial correlation')
        set(gca,'XTickLabel',[]);
        % xlabel('Time (ms)')

        %% plot BFs
        for cl = 1:2
            sp=subplot(12,3,(cl-1)*3+spl(a)+6);
            sp.Position(2) = sp.Position(2)-.004*cl;

            set(sp,'FontSize',18);
            hold on

            col = co3(cl,:);

            bf = stats.(clust{cl}).(objconds{o}).(modelnames{a}).bf;

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
            ylim(10.^(5*[-1 1]))
            sp.YTick = 10.^([-5 0 5]);
            sp.XTick = 0:100:800;
            sp.XTickLabel = [];
            xlim([min(timevect) max(timevect)+1])
            if cl == 2
                ylabel('BF')
            end
        end


        %% plot different BF
        sp=subplot(12,3,6+spl(a)+6);
        sp.Position(2) = sp.Position(2)-.004*3;
        % sp.Position(2) = sp.Position(2)+.04;
        % sp.Position(4) = sp.Position(4)-.04;
        set(sp,'FontSize',18);
        hold on

        col = [0 0 0];

        bf = stats.diff_partialcorr.(objconds{o}).(modelnames{a}).bf;

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
        ylim(10.^(5*[-1 1]))
        sp.YTick = 10.^([-5 0 5]);
        sp.XTick = 0:200:800;
        % sp.XTickLabel = [];
        xlim([min(timevect) max(timevect)+1])
        % ylabel('BF')
        xlabel('Time (ms)')

    end

    fn = sprintf('figures/FigureS%d_imageRSA_%s_only',o+2,objconds{o});
    set(gcf, 'Renderer', 'painters')
    print(gcf, '-dpng', '-r500', fn)
    im = imread([fn '.png']);
    [i,j] = find(mean(im,3) < 255); margin = 2;
    imwrite(im(min(i)-margin:max(i)+margin, min(j)-margin:max(j)+margin, :), [fn '.png'], 'png')


end