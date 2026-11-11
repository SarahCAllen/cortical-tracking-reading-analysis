%% script to clean and prepare behavioural data

clear all

% setup script here
xCOMP = 'Sarah';
addToolbox = 2; % 1 is just fieldtrip, 2 includes mi function

setup_script;

%% load and preprocess behavioural data

raw_dirs = {'bheavioural_data.csv'};

subid=[];
RT_sel = [];


for isample = 1:length(raw_dirs)

    % read in raw data
    read_time_raw = readtable([DIR.ana raw_dirs{isample}]);

    % Remove rows where target is "roble" (incorrectly listed as nonword)
    read_time_raw(strcmpi(string(read_time_raw.target), "roble"), :) = [];

    % select specific columns
    RT_temp = read_time_raw(:,["ParticipantPublicID" "ParticipantPrivateID" "ParticipantStatus" "ScreenName" "ZoneType" "ReactionTime" ...
        "Response" "Correct" "Incorrect" "ANSWER" "display"]);
    RT_sel = [RT_sel; RT_temp];

    % clean up
    clear RT_temp read_time_raw isample
    
end


%% get reaction times needed for analysis (clean)

results = [];

for isub = 1:nsub
    %reaction times
    results(isub).RTwords = RT_sel.ReactionTime(RT_sel.ParticipantPublicID==goodsub(isub) ...
        & ismember(RT_sel.display,'trials') & RT_sel.Correct == 1 & strcmp(RT_sel.ANSWER,'word'));
    results(isub).RTnonwords = RT_sel.ReactionTime(RT_sel.ParticipantPublicID==goodsub(isub) ...
        & ismember(RT_sel.display,'trials') & RT_sel.Correct == 1 & strcmp(RT_sel.ANSWER,'nonword'));

    % Save the current isub and the corresponding reaction times into results
    results(isub).isub = goodsub(isub);

    % accuracy
    results(isub).accwords = sum(RT_sel.ParticipantPublicID == goodsub(isub) ...
        & ismember(RT_sel.display, 'trials') & RT_sel.Correct == 1 & strcmp(RT_sel.ANSWER,'word'));
    results(isub).accnonwords = sum(RT_sel.ParticipantPublicID == goodsub(isub) ...
        & ismember(RT_sel.display, 'trials') & RT_sel.Correct == 1 & strcmp(RT_sel.ANSWER,'nonword'));
    
    n_correct = results(isub).accwords + results(isub).accnonwords;
    n_total = 89;

    results(isub).perc_incorrect = ...
        (n_total - n_correct) / n_total * 100;


end
