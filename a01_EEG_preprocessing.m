%% preprocessing script

% >> change ALL directories
clear all

%change the following parameters
xCOMP = 'Sarah'; 
addToolbox = 1;     % 1 for fieldtrip
setup_script;


%% load the raw data & filter
% all automatic (except for manual requirement to change sub number)

%--------------------------------------------------------------------------
isub =1;   % CHANGE SUBJECT NUMBER HERE!!! FOR EXAMPLE isub = 4; isub = 10; etc
%--------------------------------------------------------------------------
     
file = sprintf('MusicPref_sub%d.bdf',isub);     % build together the filename

% read in the two conditions separately
for icond = 1:2
    
    cfg                         = [];
    cfg.datafile                = [DIR.data file];
    cfg.trialfun                = 'ft_trialfun_general';
    cfg.trialdef.eventtype      = 'STATUS';
    cfg.trialdef.eventvalue     = usedtrg{icond,1};
    cfg.trialdef.prestim        = 2;
    cfg.trialdef.poststim       = usedtrg{icond,3};           % cut out stim duration +/-2 s for filtering/downsampling
    
    
    trl = ft_definetrial(cfg);
    
    datatmp{icond} = ft_preprocessing(trl);    
end

% merge conditions into one dataset
cfg                 = [];
cfg.keepsampleinfo  = 'ifmakessense';
data = ft_appenddata(cfg, datatmp{:});
clear datatmp

% re-reference to Cz  
cfg             = [];
cfg.reref       = 'yes';
cfg.channel     = 'all';
cfg.refchannel  = {'Cz'};

[data] = ft_preprocessing(cfg, data);

% filter data
cfg             = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.2;
cfg.hpfilttype  = 'but';
cfg.hpfiltord   = 4;
cfg.hpfiltdir   = 'twopass';

cfg.lpfilter    = 'yes';
cfg.lpfreq      = 60;
cfg.lpfilttype  = 'but';
cfg.lpfiltord   = 4;
cfg.lpfiltdir   = 'twopass'; 

cfg.detrend     = 'no';
data            = ft_preprocessing(cfg, data);

%% find and interpolate noisy channels
% manual

% find electrode neighbours for Biosemi Layout
cfg_neighb = [];
cfg_neighb.method   = 'triangulation';
cfg_neighb.layout   = 'biosemi64.lay';
cfg.channel         = {'EEG'};          % channels for which neighbours should be found
cfg.feedback        = 'yes'; 
cfg_neighb = ft_prepare_neighbours(cfg_neighb, data);

% get the electrode layout
cfg = [];
cfg.layout = 'biosemi64.lay';
layout = ft_prepare_layout(cfg, data);
 
% look at data
cfg             = [];
cfg.method      = 'summary';  
cfg.channel     = {'EEG'};
cfg.latency     = 'all';
cfg.keepchannel = 'repair';       
cfg.layout      = layout;
cfg.trials      = 'all';
cfg.neighbours  = cfg_neighb;
cfg.alim        = 1e-12;
data            = ft_rejectvisual(cfg,data);

%% re-reference to average

cfg             = [];
cfg.reref       = 'yes';
cfg.channel     = {'EEG'};
cfg.refchannel  = {'EEG'}; 
data = ft_preprocessing(cfg, data);


%% ICA 
cfg = [];
cfg.method      = 'runica';
cfg.runica.pca  = 30;
cfg.channel     = {'all'};
cfg.trials      = 'all';
comp = ft_componentanalysis(cfg,data);
ncomp = length(comp.label);

%% visualise components (involves 2 steps)

% 1) look at topography 
cfg = [];
cfg.component = 1:ncomp;           % specify the component(s) that should be plotted
cfg.layout    = 'biosemi64.lay';    % specify the layout file that should be used for plotting
cfg.comment   = 'no';
figure('units','normalized','outerposition',[0 0 1 1])  % full-window figure
ft_topoplotIC(cfg, comp)

% save an image of the component topographies
figname = sprintf('%sICAtopo_sub%d.jpg',DIR.prepro,isub);
export_fig(figname,'-r150','-a4'); 

% >> make a note of the component numbers that look like artifacts <<

%% 2) look at components again and mark the ones you want to remove
% click on green button, becomes red after you click it
% also save the to-be-removed-components as pdf (click on respective button)

cfg = [];
cfg.layout      = 'biosemi64.lay';
cfg.path        = DIR.fig;
cfg.prefix      = sprintf('ICA_sub%d',isub);    % prefix of the pdf files
cfg.powscale    = 'lin';
cfg.outputfile  = sprintf('%sICAcomps_sub%d',DIR.prepro,isub); % MAT file which contains indices of all components to reject

rejcomp = ft_icabrowser(cfg,comp);

%% REMOVE MARKED COMPONENTS 
% no need to change anything because you've marked the to-be-removed
% components in the GUI

cfg = [];
cfg.component = find(rejcomp);              % to be removed component(s)
data = ft_rejectcomponent(cfg, comp, data);


%% look at topography (you can do this at any stage of the preprocessing)

freqlim = {[1:.5:8]; [10:.5:13]};
figure('units','normalized','outerposition',[0 0 1 1])  % full-window figure
for ifreq=1:2
    % FFT
    cfg=[];
    cfg.method      = 'mtmfft';  
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.keeptrials  = 'no';
    cfg.foi         = freqlim{ifreq};    
    freq = ft_freqanalysis(cfg,data);
        
    % plot
    cfg=[];
    cfg.parameter       = 'powspctrm';
    cfg.channel         = 'all';
    cfg.layout          = 'biosemi64.lay';
    cfg.marker          = 'labels';
    cfg.zlim            = 'zeromax';
    cfg.style           = 'straight';
    %cfg.interactive     = 'on';
    subplot(1,2,ifreq)
    ft_topoplotTFR(cfg,freq);
end

% save an image of the clean topography
figname = sprintf('%stopo_clean_sub%d.jpg',DIR.prepro,isub);
export_fig(figname,'-r150','-a4');  

%% save individual conditions
% saved under cond1 (music) cond2 (story)

data2 = data;

for icond = 1:2
    
    cfg = [];
    cfg.trials = data2.trialinfo==usedtrg{icond,1};
    data = ft_selectdata(cfg,data2);
    
    save(sprintf('%sdata_clean_cond%d_sub%d',DIR.clean,icond,isub),'data','usedtrg');
end

clear data data2
    
    