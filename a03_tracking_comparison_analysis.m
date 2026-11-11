%% script for comparing speech and music tracking


clear all

% setup script here
xCOMP = 'Sarah'; 
addToolbox = 2; % 1 is just fieldtrip, 2 includes mi function

setup_script; % this loads directories, file paths and toolboxes


%% load MI data 
nfreq = 64; 
nperm = 50;
MIall = NaN(nchan,nfreq,ncond,nsub);
for icond = 1:2 
    for isub = 1:32
        load(sprintf('%scond%d_sub%d_MI_%dperm', DIR.midata, icond, isub, nperm),'realMIz','fw','nfreq');
        MIall(:,:,icond,isub) = realMIz;
    end
end

%% cluster permutation script adapted from Chalas et al.

load('cfg_neighb','cfg_neighb');

for ilabel = 1:length(cfg_neighb)

    labeltmp{ilabel} = cfg_neighb(ilabel).label;

end

freq=[];
freq.label=labeltmp;
freq.freq=fw;
freq.powspctrm=[];
freq.dimord='chan_freq';
goodsubj=1:32;
nsubj=length(goodsubj);
dat1={}; 
dat2={}; 


 counter=1;
 for ss=goodsubj
     dat1{counter}=freq;
     dat1{counter}.powspctrm=squeeze(MIall(:,:,1,ss)); % MI icond 1
     dat1{counter}.label=labeltmp';
     dat2{counter}=freq;
     dat2{counter}.label=labeltmp';
     dat2{counter}.powspctrm=squeeze(MIall(:,:,2,ss)); % MI icond 2
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
% required for a selected sample to be included % in the clustering algorithm (default=0)
cfg.tail = 0; 
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


%% topographies of average 
load('fakefreq.mat'); % load dummy fieldtrip frequency structure
load('cfg_neighb','cfg_neighb') % load neighbourhood structure

clustlist = unique(stat.posclusterslabelmat);

for iclust = clustlist(clustlist > 0)'
     
%   check if this cluster is significant
    tmp = stat.posclusterslabelmat.*stat.mask; 

%     % check how many frequency bins the cluster covers
    freqbins = sum(sum(tmp == iclust) > 0); 
    fprintf('\nCluster #%d covers %d frequency bins\n', iclust, freqbins);

    freqrange = stat.freq(sum(tmp == iclust) > 0);
    fprintf('Frequency range: %0.2f - %0.2f\n', freqrange(1), freqrange(end));

    tmax = max(abs(stat.stat(tmp == iclust)));
    d = 2*tmax / sqrt(nsub);

    % -------------------------------------------------------------------------
    % only plot the ones that cover >= 3 freq bins
    if any(any(tmp==iclust)) && freqbins >=3 

        avgX = mean(stat.stat(:,logical(sum(stat.posclusterslabelmat == iclust))),2);

        freq = [];
        freq.label = {cfg_neighb.label};
        freq.freq = 1;  % dummy freq
        freq.powspctrm = avgX;  % average of significant tstat
        freq.dimord = 'chan_freq';

        cfgp = [];
        cfgp.parameter          = 'powspctrm';
        cfgp.channel            = 'all';
        cfgp.layout             = 'biosemi64.lay';
        cfgp.zlim                = [-5 5];
        cfgp.colormap           = flipud(brewermap(128, 'RdBu'));
        cfgp.style              = 'both';
        cfgp.marker             = 'off';
        cfgp.comment            = 'no';
        cfgp.highlight          = 'on';
        cfgp.highlightchannel   = find(sum(tmp == iclust,2));
        cfgp.highlightsymbol    = 'o';
        cfgp.highlightsize      = 13;
        cfgp.highlightcolor     = [1 1 1];

        figure;
        ft_topoplotTFR(cfgp, freq);
            
       set(gcf,'renderer','painters')
       
    end
end


