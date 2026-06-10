
%% plot stimuli

ims = dir('~/Dropbox/PROJECTS/hemiraptor/hemifield-eeg/stimuli/stimuli/stim*.png');

labels = {'fish' 'bird' 'face' 'boats' 'tree' 'tools'};

words={};
for l= 1:length(labels)
    words{(l-1)*2+1} = [labels{l} '_lowercase.png'];
    words{(l-1)*2+2} = [labels{l} '_uppercase.png'];
end
    
stims = [ims(:).name words];


%% plot

f=figure(1);clf;
f.Position = [f.Position(1:2) 1500 600];
f.Resize = 'off';
aw = 70; % image size
left=50;
bufferw=2;
bufferh=2;
top = f.Position(4)-80;

coords = [1 101 201];

x=0;xy=[];
for cats = 1:6 % rows
    for i = 1:4 % columns
        x=x+1;
        if ismember(x,[101 201])
            top = top - 50;
        end
        I = imread(sprintf('%s/%s',ims(1).folder,stims{x}));
        imbox = [left+(i-1)*(aw+bufferw) top-cats*(aw+bufferh) aw aw];
        if ismember(x,coords)
            xy(end+1,:) = imbox(1:2); % mark position for later bounding box
        end
        a = axes('Units','pixels','Position',imbox,'Visible','off');
        a.XLim = [.5 250+.5];
        a.YLim = a.XLim;

        h = imshow(I);
    end
end

%
x=24;xy=[];
for cats = 1:6 % rows
    for i = 1:2 % cols
        x=x+1;

        I = imread(sprintf('%s/%s',ims(1).folder,stims{x}));
        if size(I,3)==1
            I = cat(3,I,I,I);
        end

        imbox = [(aw+bufferw)*4+left+(i-1)*(aw+bufferw) top-cats*(aw+bufferh) aw aw];

        a = axes('Units','pixels','Position',imbox,'Visible','off');
        a.XLim = [.5 250+.5];
        a.YLim = a.XLim;

        h = imshow(I(:,:,1:3));
    end
end


%% save
fn = 'Fig_stimuli';
saveas(gcf,fn,'png')
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');
