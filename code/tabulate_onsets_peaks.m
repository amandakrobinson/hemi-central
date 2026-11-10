%% tabulate onsets and peaks
load('./results/stats_decoding_central.mat')

hemi = {'L_elecs' 'R_elecs'};
heminames = {'Left' 'Right'};

%% get onsets and peaks for each hemi 

allinfo=struct();
x=0;
for h = 1:length(hemi)

    s = stats.(hemi{h});

    x=x+1;
    allinfo.Hemisphere{x,1} = heminames{h};
    allinfo.Onset(x,1) = s.onset;
    allinfo.Onset_lowCI(x,1) = s.onsetci(1);
    allinfo.Onset_highCI(x,1) = s.onsetci(2);
    allinfo.Peak(x,1) = s.peak;
    allinfo.Peak_lowCI(x,1) = s.peakci(1);
    allinfo.Peak_highCI(x,1) = s.peakci(2);
    allinfo.PeakAcc(x,1) = max(s.mu);

end

allinfo = struct2table(allinfo);

save('results/onsets_peaks.mat','allinfo');


%% tabulate onsets and peaks for image RSA
load('results/stats_imagecorrs.mat')

% get onsets and peaks for each feature x hemi

contrast = {'partialcorr_L' 'partialcorr_r'};

rsainfo=struct();
x=0;
for h = 1:length(contrast)
    for m = 1:length(modelnames)
        s = stats.(contrast{h}).(modelnames{m});

        x=x+1;
        rsainfo.Hemisphere{x,1} = heminames{h};
        rsainfo.Feature{x,1} = modelnames{m};
        rsainfo.Onset(x,1) = s.onset;
        rsainfo.Onset_lowCI(x,1) = s.onsetci(1);
        rsainfo.Onset_highCI(x,1) = s.onsetci(2);
        rsainfo.Peak(x,1) = s.peak;
        rsainfo.Peak_lowCI(x,1) = s.peakci(1);
        rsainfo.Peak_highCI(x,1) = s.peakci(2);
        rsainfo.PeakCorr(x,1) = max(s.mu);

    end
end

rsainfo = struct2table(rsainfo);

save('results/onsets_peaks.mat','allinfo','rsainfo');