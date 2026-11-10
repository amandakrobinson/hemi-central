%% run behavioural analysis

acc=[];rt=[];
for partid = 1:20

    eventsfncsv = sprintf('../hemifield-eeg/data/sub-%02i/eeg/sub-%02i_task-rsvp_events.csv',partid,partid);

    opts = detectImportOptions(eventsfncsv);
    opts = setvartype(opts, 'istarget', 'char'); 
    eventlist = readtable(eventsfncsv, opts);

    r = find(ismember(eventlist.istarget,{'L','R','C'}))';

    colldat=[];
    for a = 1:length(r) % for each target
        dat = eventlist(r(a):(r(a)+10),:); % look for responses within the next 10 stimuli (2 seconds)

        cond = dat.condition(1);
        targpos = dat.istarget(1);
       
        if contains(targpos,'L')
            tpos = 1;
        elseif contains(targpos,'C')
            tpos=2;
        elseif contains(targpos,'R')
            tpos=3;
        end

        if dat.stimL(1)>-1
            spos = 1;
        elseif dat.stimC(1)>-1
            spos = 2;
        elseif dat.stimR(1)>-1
            spos = 3;
        end

        idx = find(dat.response==1);
        if length(idx)==1
            ridx = idx;
        elseif length(idx)>1
            ridx = idx(1);
        else
            ridx = 1;
        end

        rdat = dat(ridx,:);

        colldat(a,:) = [cond tpos spos rdat.response rdat.rt ridx];

    end

    % now collate accuracy and RT per condition x position of target
    for cond = 1:3
        % results per target position
        for pos = 1:3
            id = find(colldat(:,1)==cond&colldat(:,2)==pos);

            acc(cond,pos,partid) = sum(colldat(id,4))/length(id);
            rt(cond,pos,partid) = median(colldat(id,5));
        end
        idx = find(colldat(:,1)==cond);
        meancond(cond,partid) = sum(colldat(idx,4))/length(idx);
        ntargscond(cond,partid) = length(idx);


    end
    ntargs(partid) = length(colldat);
    allacc(partid) = sum(colldat(:,4))/length(colldat);
    allrt(partid) = median(colldat(colldat(:,4)==1,5));

end

%% plot results

mu = mean(meancond,2);
se = std(meancond,[],2)'/sqrt(size(meancond,2));

cols = plasma(3);
cols(1,:) = [0 0 0];

figure
set(gcf,'Position',[1838 155 1547 593])
subplot(1,2,1)
set(gca,'FontSize',20)
hold on
ylim([0.5 1])
xlim([.5 3.5])
for c = 1:size(mu,1) % for each cond
    bar(c,mu(c),'FaceColor',cols(c,:))
    errorbar(c,mu(c),se(c),'Color',cols(c,:),'HandleVisibility','off')
end
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'Central' 'Single peripheral','Dual peripheral'})
xlabel('Stimulus condition')
ylabel('Proportion correct')
title('Overall results')

fn = 'figures/FigureS1_behaviouralaccuracy';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

