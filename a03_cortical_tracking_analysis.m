%% script for analysing cortical tracking of speech and music

clear all

% setup script here
xCOMP = 'Sarah'; 
addToolbox = 2; % 1 is just fieldtrip, 2 includes mi function

setup_script; % this loads directories, file paths and toolboxes

%% load MI (or also conditional MI)

nfreq = 64; 
nperm = 50;
MIall = NaN(nchan,nfreq,ncond,nsub);
for icond = 1:2 
    for isub = 1:32
        load(sprintf('%scond%d_sub%d_MI_%dperm', DIR.midata, icond, isub, nperm),'realMIz','fw','nfreq');
        MIall(:,:,icond,isub) = realMIz;
    end
end


%% Nikos cluster permutation script with some changes

% cluster permutation stats
load('cfg_neighb','cfg_neighb');


for ilabel = 1:length(cfg_neighb)

    labeltmp{ilabel} = cfg_neighb(ilabel).label;

end

icond = 2; %1 for music, 2 for speech

freq=[];
freq.label=labeltmp;
freq.freq=fw;
freq.powspctrm=[];
freq.dimord='chan_freq';
goodsubj=1:32;
nsubj=length(goodsubj);
dat1={}; %  MI per channel
dat2={}; % zeros


 counter=1;
 for ss=goodsubj
     dat1{counter}=freq;
     dat1{counter}.powspctrm=squeeze(MIall(:,:,icond,ss));
     dat1{counter}.label=labeltmp';
     dat2{counter}=freq;
     dat2{counter}.label=labeltmp';
     dat2{counter}.powspctrm = zeros(size(dat1{counter}.powspctrm));
     counter=counter+1;
 end

 %% do analysis

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster'; %cluster
cfg.parameter = 'powspctrm';
cfg.alpha = .05;
cfg.clusteralpha = .05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan = 3; % minimum number of neighborhood freq channels that is
% required for a selected sample to be included in the clustering algorithm (default=0)
cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.neighbours=[];
cfg.neighbours=cfg_neighb;
cfg.numrandomization = 1000;
% set up design matrix
cfg.design = [1:nsubj 1:nsubj; ones(1, nsubj) 2*ones(1, nsubj)];
cfg.uvar = 1;
cfg.ivar = 2;
% run t-test
stat = ft_freqstatistics(cfg, dat1{:}, dat2{:});

tstat=stat.mask.*stat.stat;


%% topographies of average (significant cluster)

% load dummy fieldtrip frequency structure
load('fakefreq.mat');
% load neighbourhood structure
load('cfg_neighb','cfg_neighb')

% choose positive / negative clusters to plot
clustlist = unique(stat.posclusterslabelmat);
% clustlist = unique(stat.negclusterslabelmat);


for iclust = clustlist(clustlist > 0)'

%   check if this cluster is significant
    tmp = stat.posclusterslabelmat.*stat.mask;    
   

%    check how many frequency bins the cluster covers
    freqbins = sum(sum(tmp == iclust) > 0); 
    freqrange = stat.freq(sum(tmp == iclust) > 0);
    fprintf('\nCluster #%d covers %d frequency bins\n', iclust, freqbins);
    fprintf('Frequency range: %0.2f - %0.2f\n', freqrange(1), freqrange(end));

    mask = (tmp == iclust);           % logical mask (channels×freqs)

    % get absolute t-values inside cluster and linear index of peak
    tvals = abs(stat.stat(mask));
    [tmax, linidx] = max(tvals);

    % get frequency indices corresponding to true mask entries
    [~, freq_idx] = find(mask);

    peak_freq = stat.freq(freq_idx(linidx));

    fprintf('tmax = %g at frequency = %g Hz\n', tmax, peak_freq);

    tmax = max(abs(stat.stat(tmp == iclust)));
    d = 2*tmax / sqrt(nsub);
    fprintf('tmax = %0.2f, Cohen d = %0.2f\n', tmax, d);

    % -------------------------------------------------------------------------
    % only plot the ones that cover >= 3 freq bins
    if any(any(tmp==iclust)) && freqbins >=3 

        avgX = mean(stat.stat(:,logical(sum(stat.posclusterslabelmat == iclust))),2);

        freq = [];
        freq.label = {cfg_neighb.label};
        freq.freq = 1;  % dummy freq, necessary to plot
        freq.powspctrm = avgX;  % average of significant tstat
        freq.dimord = 'chan_freq';

        cfgp = [];
        cfgp.parameter          = 'powspctrm';
        cfgp.channel            = 'all';
        cfgp.layout             = 'biosemi64.lay';
        cfgp.zlim                = [0 5];
        cfgp.colormap           = crameri('-devon'); 
        cfgp.style              = 'both';
        cfgp.marker             = 'off';
        cfgp.comment            = 'no';
        cfgp.highlight          = 'on';
        cfgp.highlightchannel   = find(sum(tmp == iclust,2));
        cfgp.highlightsymbol    = 'o';
        cfgp.highlightsize      = 5;

        figure;
        ft_topoplotTFR(cfgp, freq);
        colorbar;

         % then change the marker fills to white
        h = findobj(gcf,'Marker','o');          
        set(h,'MarkerFaceColor','w','MarkerEdgeColor','black');


    end
end

