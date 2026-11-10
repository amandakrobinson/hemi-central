% Group-level signal quality analysis for the LH and RH clusters
% Run run_snr(subjectnr) for every participant first.
%
% Outputs:
%   results/snr_all.csv     - long table, one row per subject x cluster x window
%   results/snr_summary.csv - RMS(SME), mean SNR and LH-RH comparisons
%   figures/snr_erps.png    - group ERPs per cluster, RMS(SME)
%
% Stats: paired LH vs RH comparisons using bayes factors

subjects = 1:20;

%% load
T = table();
for s = subjects
    d = load(sprintf('results/sub-%02i_snr.mat', s));
    T = [T; d.snr]; %#ok<AGROW>
    for c = {'LH','RH'}
            E.(c{1})(s,:) = d.erp.(c{1});
    end
end
timevect = d.timevect;
windows  = d.windows;
writetable(T, 'results/snr_all.csv');

%% summary: RMS(SME) per cluster, and paired LH vs RH tests
summ = {};

for w = unique(T.window)'
    L = T(strcmp(T.cluster,'LH') & strcmp(T.window,w), :);
    R = T(strcmp(T.cluster,'RH') & strcmp(T.window,w), :);
    L = sortrows(L,'subject'); R = sortrows(R,'subject');
    assert(isequal(L.subject, R.subject))

    rmsL = sqrt(mean(L.aSME.^2));   rmsR = sqrt(mean(R.aSME.^2));

    stat_sme = bayesfactor_R_wrapper((L.aSME-R.aSME)',...
            'returnindex',1,'verbose',false,...
            'args','mu=0,rscale="medium"');   % BF10, paired
    stat_snr = bayesfactor_R_wrapper((L.SNR-R.SNR)',...
            'returnindex',1,'verbose',false,...
            'args','mu=0,rscale="medium"');   % BF10, paired;

    summ(end+1,:) = {w{1}, rmsL, rmsR, mean(L.SNR), mean(R.SNR), ...
        stat_sme, stat_snr}; %#ok<SAGROW>
end

S = cell2table(summ, 'VariableNames', {'window', ...
    'RMS_aSME_LH','RMS_aSME_RH','meanSNR_LH','meanSNR_RH', ...
    'bf_aSME', 'bf_SNR'});
writetable(S, 'results/snr_summary.csv');
disp(S)

%% figure: group ERPs per cluster and stimulus type
cols = [0.40 0.15 0.45; 0.10 0.60 0.55]; % LH purple, RH teal (as in paper)
figure('Position',[77 211 1436 464]);
subplot(1, 3, 1); hold on
for c = 1:2
    cn = {'LH','RH'};
    x  = E.(cn{c})(subjects,:);
    m  = mean(x,1);
    se = std(x,[],1) / sqrt(size(x,1));
    fill([timevect fliplr(timevect)], [m-se fliplr(m+se)], cols(c,:), ...
        'FaceAlpha', .2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(timevect, m, 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', cn{c});
end
ylim([-20 20])

yl = ylim;
for w = 1:numel(windows)
    patch(windows(w).range([1 2 2 1]), yl([1 1 2 2]), [.5 .5 .5], ...
        'FaceAlpha', .08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(windows(w).range(1)+10,-17,windows(w).name)
end

xline(0, 'HandleVisibility', 'off');
xlim([-100 500]); title('Mean ERPs'); xlabel('Time (ms)');
ylabel('Amplitude (CSD units)'); legend();
set(gca,'FontSize',22)

subplot(1,3,2); hold on
cn   = {'LH','RH'};
wn   = {windows.name};
nw   = numel(wn);
ns   = numel(subjects);
off  = [-.15 .15];                 % x offset for LH / RH within a window

% subject x window x cluster, averaged over stimulus types
sme = nan(ns, nw, 2);
for c = 1:2
    for w = 1:nw
        for i = 1:ns
            idx = strcmp(T.cluster, cn{c}) & strcmp(T.window, wn{w}) & ...
                  T.subject == subjects(i);
            sme(i,w,c) = T.aSME(idx);
        end
    end
end

% paired subject lines
for w = 1:nw
    plot(w + off, squeeze(sme(:,w,:))', '-', 'Color', [.8 .8 .8], ...
        'HandleVisibility', 'off');
end

% individual points + group mean ± SE
for c = 1:2
    x = repmat((1:nw) + off(c), ns, 1);
    y = sme(:,:,c);
    scatter(x(:), y(:), 25, cols(c,:), 'filled', 'MarkerFaceAlpha', .35, ...
        'HandleVisibility', 'off');
    m  = sqrt(mean(y.^2, 1, 'omitnan'));
    ci = nan(2, nw);
    for w = 1:nw
        yw = y(~isnan(y(:,w)), w);
        ci(:,w) = bootci(2000, @(v) sqrt(mean(v.^2)), yw);
    end
    errorbar((1:nw) + off(c), m, m - ci(1,:), ci(2,:) - m, 'o', ...
        'Color', cols(c,:), 'MarkerFaceColor', cols(c,:), 'MarkerSize', 9, ...
        'LineWidth', 2, 'CapSize', 0, 'DisplayName', cn{c});
end

% BF10 (LH vs RH aSME) above each window
yl = ylim;
for w = 1:nw
    bf = S.bf_aSME(strcmp(S.window, wn{w}));
    text(w, yl(2), sprintf('BF = %.3g', bf), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 14);
end
ylim([0 10])
xlim([.5 nw+.5]); xticks(1:nw); xticklabels(wn);
xlabel('Window'); ylabel('aSME (CSD units)'); 
title('RMS(aSME) per window');
legend('Location', 'best');
set(gca, 'FontSize', 22)

%
subplot(1,3,3); hold on

% subject x window x cluster, averaged over stimulus types
snr = nan(ns, nw, 2);
for c = 1:2
    for w = 1:nw
        for i = 1:ns
            idx = strcmp(T.cluster, cn{c}) & strcmp(T.window, wn{w}) & ...
                  T.subject == subjects(i);
            snr(i,w,c) = T.SNR(idx);
        end
    end
end

% paired subject lines
for w = 1:nw
    plot(w + off, squeeze(snr(:,w,:))', '-', 'Color', [.8 .8 .8], ...
        'HandleVisibility', 'off');
end

% individual points + group mean ± SE
for c = 1:2
    x = repmat((1:nw) + off(c), ns, 1);
    y = snr(:,:,c);
    scatter(x(:), y(:), 25, cols(c,:), 'filled', 'MarkerFaceAlpha', .35, ...
        'HandleVisibility', 'off');
    m  = mean(y, 1, 'omitnan');
    se = std(y,[],1)/sqrt(size(y,1));
    errorbar((1:nw) + off(c), m, se, 'o', ...
        'Color', cols(c,:), 'MarkerFaceColor', cols(c,:), 'MarkerSize', 9, ...
        'LineWidth', 2, 'CapSize', 0, 'DisplayName', cn{c});
end

% BF10 (LH vs RH aSME) above each window
yl = ylim;
for w = 1:nw
    bf = S.bf_SNR(strcmp(S.window, wn{w}));
    text(w, yl(2)-3, sprintf('BF = %.3g', bf), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 14);
end
% ylim([0 10])
xlim([.5 nw+.5]); xticks(1:nw); xticklabels(wn);
xlabel('Window'); ylabel('SNR'); 
title('SNR per window');
legend('Location', 'best');
set(gca, 'FontSize', 22)

%
fn = './figures/FigureS3_snr_erps';
set(gcf, 'Renderer', 'painters')
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=0;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

