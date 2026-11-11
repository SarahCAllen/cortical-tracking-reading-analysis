%% script for analysing tracking and correlations with behaviour

clear all

% setup script here
xCOMP = 'Sarah'; 
addToolbox = 2; % 1 is just fieldtrip, 2 includes mi function

setup_script; % this loads directories, file paths and toolboxes

%% load MI

nfreq = 64; 
nperm = 50;
MIall = NaN(nchan,nfreq,ncond,nsub);
for icond = 1:2 
    for isub = 1:32
        load(sprintf('%scond%d_sub%d_MI_%dperm', DIR.midata, icond, isub, nperm),'realMIz','fw','nfreq');
        MIall(:,:,icond,isub) = realMIz;
    end
end

%% load behavioural data (already pre-processed in script a01_clean_behavioural_data.m)

load("clean_behavioural_task", "results", "RTs", "ACCs", "RTmedian", "ntrials");


%% check if any participants have too many trials excluded
% ntrials: nsub x ncond numeric matrix

% combine conditions per subject
combined = sum(ntrials, 2); 

% stats and z-scores
mu_comb    = mean(combined);
sigma_comb = std(combined, 0);
z_comb     = (combined - mu_comb) / sigma_comb;

% participants outside ±2 SD
outside_subj = find(abs(z_comb) > 2); 

% print results
fprintf('Combined mean=%.2f, sd=%.2f\n', mu_comb, sigma_comb);
fprintf('Number of participants outside ±2 SD: %d\n', numel(outside_subj));
if ~isempty(outside_subj)
    fprintf('Participant numbers: %s\n', mat2str(outside_subj.'));
end

%% removing RT outliers for correlation below

keep = setdiff(1:nsub, outside_subj);
RTs =RTs(keep,:);
RTmedian =RTmedian(keep,:);
MIall = MIall(:,:,:,keep);
nsub = numel(keep);

%% read in dprime (for Supplementary)

dprime = readtable("dprime.csv");

%% 01) script for correlation (cluster-based permutation - custom-made)

corrtype = 'Pearson'; 

% 1 = music  2 = speech
icond = 1; 

ibeh = 1; 

if ibeh == 1
    beh_meas = mean(RTmedian,2);    % reaction times median overall (main results)
elseif ibeh == 2
    beh_meas = RTmedian(:,1);       % reaction times median words (supplementary)
elseif ibeh == 3
    beh_meas = RTmedian(:,2);       % reaction times median nonwords (supplementary)
elseif ibeh == 4
    beh_meas = ACCs(:,1);           % accuracy words (supplementary)
elseif ibeh == 5
    beh_meas = ACCs(:,2);           % accuracy non-words (supplementary)
elseif ibeh == 6
    beh_meas = dprime(:,1);         % dprime (supplementary)
end

load('cfg_neighb','cfg_neighb')

% true rho-values (correlation for all chans individually)
clear RHO PVAL 
RHO = NaN(nchan,1);
for ichan = 1:nchan
    for ifreq = 1:nfreq
    [RHO(ichan,ifreq),PVAL(ichan,ifreq)]=corr(squeeze(MIall(ichan,ifreq,icond,:)),beh_meas, 'Type', corrtype);
    end
end

X = corr2z(RHO);

% make random b-matrix by shuffling participant labels
nperm = 1000; 
XR = NaN(nchan,nfreq, nperm);
for iboot=1:nperm
    psub = randsample(nsub,nsub);
    
    RHOp = NaN(nchan,nfreq);
    for ichan = 1:nchan          
        for ifreq = 1:nfreq
        [RHOp(ichan,ifreq),PVAL(ichan,ifreq)]=corr(squeeze(MIall(ichan, ifreq, icond, psub)),beh_meas, 'Type', corrtype);                
        end
    end

    XR(:,:,iboot) = corr2z(RHOp);      % z-transformed r-vals    
end


%% permutation test
cfg=[];
cfg.neighbours = cfg_neighb;

% for z-scored rvals
cfg.critvaltype = 'prctile';
q = 0.05; 

q = q*100/2;
cfg.critval = [q, 100-q];  % two-sided

cfg.conn = 8; % connectivity criterion
cfg.clusterstatistic = 'maxsum';
cfg.minsize = 3; % minimal cluster size 
cfg.minfreq = 3; % frequency threshold for cluster
cfg.pval = 0.05; % threshold to select significant clusters
cfg.df = nsub-1;

[PosClus,NegClus] = eegck_clusterstats_eeg(cfg,X,XR);

cfg_stat = [];
cfg_stat.p = .05;
cfg_stat.fw = fw;
AK_clustinfo(PosClus,NegClus,cfg_stat);


%% topography of averaged X

load('fakefreq.mat'); % load dummy fieldtrip frequency structure
load('cfg_neighb','cfg_neighb') % load neighbourhood structure

plotclust = [1];

% choose which cluster to check
for iclust = plotclust


       % find max r and use for Cohen's d analysis
        rmax = max(abs(X(NegClus.maskSig==iclust)));
        d = 2*rmax / sqrt(1-rmax^2);
  
        avgX = mean(X(:,logical(sum(xeffect.maskSig == iclust))),2);

        freq = [];
        freq.label = {cfg_neighb.label};
        freq.freq = 1;  % dummy freq
        freq.powspctrm = avgX;  % average of significant X
        freq.dimord = 'chan_freq';

        cfgp = [];
        cfgp.parameter          = 'powspctrm';
        cfgp.channel            = 'all';
        cfgp.layout             = 'biosemi64.lay';
        cfgp.zlim               = 'maxabs';
        cfgp.colormap           = brewermap(128, '*RdBu');
        cfgp.style              = 'both';
        cfgp.marker             = 'off';
        cfgp.comment            = 'no';
        cfgp.highlight          = 'on';
        cfgp.highlightchannel   = find(sum(xeffect.maskSig == iclust,2));
        cfgp.highlightsymbol    = 'o';
        cfgp.highlightsize      = 8;
        cfgp.highlightcolor     = [1 1 1];

        figure;
        ft_topoplotTFR(cfgp, freq);
        colorbar;

end







