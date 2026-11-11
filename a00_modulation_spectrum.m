%% load folders
clear all
% >> add setup script
xCOMP = 'Sarah'
setup_script;

%% Modulation spectrum
Fs = soundstim.Fs

% load soundstim here
figure
colors = lines(2);
for icond=1:2
    

y = soundstim.soundstim{1,icond};

nsec = 20;  
segl = nsec*Fs; 
nseg = ceil(length(y)/segl);

% zero-pad signal at the end
lng = length(y);
ytmp = zeros(nseg*segl,1);
ytmp(1:lng) = y(:,1);   


clear ms
k=1;
for iseg = 1:nseg
    yseg = ytmp(k:k+segl-1,1);  
    % resample the recordign to 16 k Hz
    x=resample(yseg,16e3,Fs);
    % generate the auditory spectrogram
    [v,cf]=wav2aud2(x,[5 8 -2 0]); 
    % compute the narrow-band modulation spectrum
    [ms(:,iseg),f]=Modulation_Spectrum(v,200,'log');
    % increase counter
    k=k+segl;
end

% calculate the root-mean-square of the modulation spectrum
% across all the 6-second duration chunks
ms_rms=sqrt(mean(ms.^2,2));
ms_rms=ms_rms/max(ms_rms([f<32 & f>=.5]));



% Find modulation peak
valid_idx = find(f >= 0.5 & f <= 32);
[~, idx_peak] = max(ms_rms(valid_idx));
idx_peak = valid_idx(idx_peak);
peak_freq = f(idx_peak);
peak_val = ms_rms(idx_peak);


% figure; 
hold on
plot(f(1:size(ms,1)),ms_rms,'LineWidth',2, 'Color', colors(icond, :)) 


    % Mark the peak
    plot(peak_freq, peak_val, 'o', 'Color', colors(icond, :), 'MarkerFaceColor', colors(icond, :))
    text(peak_freq, peak_val * 1.05, sprintf('%.2f Hz', peak_freq), ...
        'Color', colors(icond, :), 'FontSize', 10, 'HorizontalAlignment', 'center');

    
    fprintf('Condition %d → Peak modulation frequency = %.2f Hz\n', icond, peak_freq);
end

xlim([0 20]);
xlabel('frequency (Hz)')
ylabel('normalised amplitude')
title(sprintf('Modulation spectra (%ds segments)',nsec));


