clear all

%change the following parameters
xCOMP = 'Sarah'; 
addToolbox = 2;     % 1 for fieldtrip
setup_script; % loads all paths and directories
nsurr=50;

%% start loops

    % Load previously computed lags per participant
    load('max_lags_by_subject.mat', 'max_delay_samples_all');

for icond = 1:ncond
    
    % Access relevant lags for this condition and subject
    max_delay_samples = max_delay_samples_all(:, icond);

    
    % load stimulus envelope
    load(sprintf('%s%s',DIR.ana,condstr{3,icond}),'ENV','Rate');
    nEnv = length(ENV);
    
    dEnv=[0 diff(ENV)]; %computes the first derivative of the envelope of an acoustic signal and prepends a zero to the result
    fb = cwtfilterbank('SignalLength',numel(ENV),'SamplingFrequency',150,'FrequencyLimits',[0.25 20],'WaveletParameters',[3,10]); %creates a continuous wavelet transform (CWT) filter bank fb.
    
    % frequency transform of Envelope, using filterbank
    [Envm,fw]=wt(fb,ENV);
    Envm=Envm./sum(abs(Envm),2);
    
    Envm(isnan(Envm))=0;
    % frequency transform of Envelope Difference (first derivative), using filterbank
    [dEnvm,~]=wt(fb,dEnv);
    dEnvm=dEnvm./sum(abs(dEnvm),2);
    dEnvm(isnan(dEnvm))=0;

    nfreq = length(fw);


    for isub = 1:nsub
        tic
        fprintf('\n*** Analysing participant #%d ***\n',isub)
        realMI=NaN(nchan,nfreq);
        realMIz=NaN(nchan,nfreq);
        surrMI=NaN(nchan,nfreq,nsurr);


        % load clean eeg data
        data_struct = load(sprintf('%sdata_clean_cond%d_sub%d', DIR.clean, icond, goodsub(isub)),'data');
        data = data_struct.data;  
        cfg = [];
        cfg.resamplefs  = Rate;
        cfg.detrend     = 'no';
        [data] = ft_resampledata(cfg, data);

        % find zero in data.time
        zerosamp = nearest(data.time{1},0);

        % cut out data from zero time point, keep the rest complete
        cfg = [];
        cfg.begsample = zerosamp;
        cfg.endsample = length(data.time{1});
        data = ft_redefinetrial(cfg,data);

        for ichan=1:nchan
            
            fprintf('\n*** Analysing participant %d (ID #%d) - channel %d ***\n',isub,goodsub(isub),ichan)
            data_current = data.trial{1};
            dataloop=data_current(ichan,:); % go through each channel
            datacut = dataloop(max_delay_samples(isub):max_delay_samples(isub)+nEnv-1);

            [datac1,~]=wt(fb,datacut(1,:));
            datac1=datac1./sum(abs(datac1),2);

            x1=copnorm(real(Envm)') ;
            x2= copnorm(imag(Envm)');
            x3=copnorm(real(dEnvm)') ;
            x4= copnorm(imag(dEnvm)');
            y1=copnorm(real(datac1)') ;
            y2= copnorm(imag(datac1)');


            for ifreq=1:nfreq
                realMI(ichan,ifreq)=mi_gg([x1(:,ifreq) x2(:,ifreq) x3(:,ifreq) x4(:,ifreq)],[y1(:,ifreq) y2(:,ifreq)]);
            end

            % surrogate MI values
            tmp=NaN([nsurr nfreq]);


            %using random shuffling of segments
            segsize = 750; %5 secs
            
            segm = NaN(ceil(length(x1) / segsize), 2);  % Pre-allocate 'segm' for performance
            segm(:,1) = 1:segsize:length(x1);
            segm(:,2) = segm(:,1)+segsize-1;
            segm(end,2) = length(x1);

             parfor isurr = 1:nsurr
                % randomisation of segments:
                randseg = randperm(length(segm));   
                x1r = [];
                x2r = [];
                x3r = [];
                x4r = [];
                for iseg = 1:length(segm)
                    x1r = [x1r; x1([segm(randseg(iseg),1):segm(randseg(iseg),2)],:)];
                    x2r = [x2r; x2([segm(randseg(iseg),1):segm(randseg(iseg),2)],:)];
                    x3r = [x3r; x3([segm(randseg(iseg),1):segm(randseg(iseg),2)],:)];
                    x4r = [x4r; x4([segm(randseg(iseg),1):segm(randseg(iseg),2)],:)];
                end
                for ifreq=1:nfreq
                    tmp(isurr,ifreq) = mi_gg([x1r(:,ifreq) x2r(:,ifreq) x3r(:,ifreq) x4r(:,ifreq)],[y1(:,ifreq) y2(:,ifreq)]);
                end

            end     % isurr


            sm=nanmean(tmp,1); 
            sd=nanstd(tmp,0,1);

            realMIz(ichan,:)=(squeeze(realMI(ichan,:))-sm)./sd;

            surrMI(ichan,:,:) = tmp';


        end     % ichan        

        savename = sprintf('%scond%d_sub%d_MI_%dperm', DIR.midata, icond, isub, nsurr);
        save(savename,'realMI','realMIz','surrMI','fw','nfreq','goodsub');
                
        toc
    end     %isub


end   %icond
