%% script for individual lags

clear all

% setup script here
xCOMP = 'Sarah';
addToolbox = 2; % 1 is just fieldtrip, 2 includes mi function

setup_script; % this loads directories, file paths and toolboxes

%%

condstr = {'music', 'story';'LinRoundtree_Fluid', 'elvesshoemaker'}

% frequency range to analyse
freqlim = {.25 20 '20Hz'};   % min, max, string for filename

% delays/lags to analyse (0 to 500 ms in 13.3 ms increments)
delay_samples = [301:2:375]; 
ndel = length(delay_samples);

MItrue = NaN(nsub,ncond,ndel);    % create empty matrix for 2 conditions and all subs and channels (only want to do this once or will overwrite previous)

%%
for icond = 1:ncond

    %%% process audio--------------------------------------------
    fname = sprintf('ENV_%s150Hz.mat',condstr{2,icond})
    load(fname,'ENV','Rate');

    % create derivative
    dEnv=[0 diff(ENV)];

    % make filter for audio
    ARG.mode = 'butter';
    ARG.Freqs = [freqlim{ifreq,1} freqlim{ifreq,2}];
    ARG.order = 3;
    ARG.rate = 150;
    FF = ck_filt_makefilters(ARG);

    zcopnorm = @(x) copnorm([real(x(:)) imag(x(:))]);

    % USE no OFFSET TO REMOVE EVOKED RESPONSE
    OFFSET = 0;
    ndur = length(ENV);
    maxl = ndur-OFFSET;

    % copnorm SPEECH
    storysig = ck_filt_applyfilter(ENV,FF{1});
    storysig = hilbert(storysig(OFFSET+[1:maxl])');
    storysig = hilbert(storysig([1:maxl])');
    storysig = zcopnorm(storysig(:));

    % build derivative
    storyder = ck_filt_applyfilter(dEnv,FF{1});
    storyder = hilbert(storyder([1:maxl])');
    storyder = zcopnorm(storyder(:));

    %%% start subject loop
    for isub = 1:nsub

        % load data
        load(sprintf('%sdata_clean_cond%d_sub%d',DIR.clean,icond,goodsub(isub)),'data','usedtrg');

        % highpass / lowpass filter
        cfg             = [];
        cfg.hpfilter    = 'yes';
        cfg.hpfreq      = freqlim{1};
        cfg.hpfilttype  = 'but';
        cfg.hpfiltord   = 3;
        cfg.hpfiltdir   = 'twopass';
        cfg.lpfilter    = 'yes';
        cfg.lpfreq      = freqlim{2};
        cfg.lpfilttype  = 'but';
        cfg.lpfiltord   = 3;
        cfg.lpfiltdir   = 'twopass';
        cfg.detrend     = 'no';
        data            = ft_preprocessing(cfg, data);

        % resample to 150Hz
        cfg             = [];
        cfg.resamplefs  = Rate;
        data            = ft_resampledata(cfg, data);

        % cut out only Cz
        cfg             = [];
        cfg.channel     = {'Cz'};
        dataCz          = ft_preprocessing(cfg, data);

        fprintf('\nInfo loop for: sub #%d - condition %d\n',isub,icond);

        for ilag = 1:ndel
            % cut out story *plus lag* & plus offset
            cfg             = [];
            cfg.begsample   = delay_samples(ilag);          % same lag for all
            cfg.endsample   = delay_samples(ilag)+maxl-1;
            data            = ft_redefinetrial(cfg, dataCz);

            fprintf('\n processing lag %d \n',delay_samples(ilag));
            % bookkeeping
            ntime=size(data.trial{1},2);

            % get data out of fieldtrip structure
            data2=data.trial{1};

            % copnorm EEG data
            tmpEEG = zcopnorm(hilbert(squeeze(data2(:,:))'));

            % MI (without randomisation)
            MItrue(isub,icond,ilag) = mi_gg([storysig storyder], tmpEEG,0,1);

            
        end     % lag loop

        fprintf('\n*** Sub #%d condition %d (%s) done ***\n\n',isub,icond,freqlim{ifreq,3});

    end     % sub loop
end     % cond loop

% save data
save(sprintf('MI_%s_allsubs_Cz',freqlim{ifreq,3}),'MItrue','freqlim','nsub','ncond','nchan','goodsub','delay_samples');


%% find lag at which frequency-averaged MI is highest

% Preallocate
max_vals_all = NaN(nsub, ncond);
max_delay_samples_all = NaN(nsub, ncond);

% Loop over conditions and subjects
for icond = 1:ncond
    for isub = 1:nsub
        [max_vals_all(isub, icond), max_idx] = max(MItrue(isub, icond, 5:end));
        max_delay_samples_all(isub, icond) = delay_samples(max_idx + 4); 
    end
end

% Subtract 300 from max MI values (2 sec leading time from preprocessing)
max_delay_samples_all = max_delay_samples_all - 300;

%%
% Save for later use
save('max_lags_by_subject.mat', 'max_delay_samples_all', 'max_vals_all', 'goodsub');