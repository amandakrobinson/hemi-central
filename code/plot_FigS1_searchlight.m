%% add paths
addpath('~/Dropbox/MATLAB/matplotlib/')
addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa/')
addpath('~/Dropbox/MATLAB/fieldtrip-20171104/');

ft_defaults;

%% collate searchlight results

res_cell = {};
for s=1:20
    x=load(sprintf('./results/sub-%02i_ch_searchlight_multiclass.mat',s));
    res = cosmo_slice(x.res,x.res.sa.cond==1);
    res_cell{s} = res;
end
res_searchlight = cosmo_stack(res_cell);


%% get layout
layout=cosmo_meeg_find_layout(res_searchlight,'label_threshold',.80);

cfg={};
ft = ft_timelockanalysis(cfg,cosmo_map2meeg(res_searchlight));

%% plot

f=figure(1);clf
f.Position = [f.Position(1:2) 1400 700];
co = tab10; %vega10;

tv = res_searchlight.a.fdim.values{2};
times = 60:20:340;

for t=1:length(times)
    clear a1
    a1=subplot(3,5,t);hold on

    % topo
    axis off
    tp = find(res_searchlight.a.fdim.values{2}==times(t));
    x = mean(res_searchlight.samples(:,res_searchlight.fa.time==tp),1);
    ft.avg = repmat(x',1,size(ft.avg,2));

    mrange = [1/36 .045];

    a1pos = a1.Position;
    a=axes('Position',[a1pos(1)-.03 a1pos(2)-.02 a1pos(3)+.06 a1pos(4)+.04]);
    hold on
    % co3 = [linspace(1,co(1,1),100);...
    %     linspace(1,co(1,2),100);...
    %     linspace(1,co(1,3),100)]';
    co3 = viridis;
    % show figure with plots for each sensor
    cfg = [];
    cfg.zlim = mrange;
    cfg.xlim = tv([tp tp+1]);
    cfg.layout = layout;
    % cfg.showscale = 'no';
    cfg.figure = a;
    cfg.interactive ='no';
    cfg.comment = 'no';
    % cfg.markersymbol = '.';
    cfg.style = 'straight';
    cfg.gridscale = 200;
    ft_info off verbose
    ft_topoplotER(cfg, ft);
    a.FontSize = 12;
    a.Colormap = co3;
    set(a.Children,'LineWidth',.5)
    title(sprintf('%d ms',tv(tp)))
    % ha(t,d)=a;
end

% save
fn = 'figures/FigS1_searchlight';
set(gcf, 'Renderer', 'painters')
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');