%% lateralisation script 

clear all

% setup script here
xCOMP = 'Sarah';
setup_script; % this loads directories, file paths and toolboxes

%% load data

load("results", "NegClus", "PosClus", "X", "XR", "fw"); 

%%
xeffect = NegClus;
nclust = max(unique(xeffect.maskSig));

for iclust = 3 % iclust 1 for music, iclust 3 for speech
    freqbins = sum(sum(xeffect.maskSig == iclust) > 0);
    freqrange = fw(sum(xeffect.maskSig == iclust) > 0);
end



%%  lateralisation

load('cfg_neighb','cfg_neighb')
labels = {cfg_neighb.label};

iclust = 3;
sigmask = xeffect.maskSig == iclust;

% disp(labels(any(sigmask,2))')


leftCh  = {'F1','F3','F5','F7','FC1','FC3','FC5','C1','C3','C5','CP1','CP3','CP5','P1','P3','P5','PO3','PO7','TP7','T7','P7','Fp1','AF3','AF7','O1'};
rightCh = {'F2','F4','F6','F8','FC2','FC4','FC6','C2','C4','C6','CP2','CP4','CP6','P2','P4','P6','PO4','PO8','TP8','T8','P8','Fp2','AF4','AF8','O2'};
midCh   = {'Fz','FCz','Cz','CPz','Pz','POz','Oz','AFz','Fpz','Iz'};

aE = []; cE = []; nLeft = 0; nRight = 0; nMid = 0;
for e = 1:numel(labels)
    if ~any(sigmask(e,:)), continue; end
    li = find(strcmp(leftCh,  labels{e}));
    ri = find(strcmp(rightCh, labels{e}));
    isMid = false;
    if     ~isempty(li), partner = rightCh{li}; nLeft  = nLeft  + 1;
    elseif ~isempty(ri), partner = leftCh{ri};  nRight = nRight + 1;
    elseif any(strcmp(midCh, labels{e})), partner = labels{e}; isMid = true; nMid = nMid + 1;
    else,  continue;                           % any other unpaired electrode
    end
    p = find(strcmp(labels, partner));
    if isempty(p), continue; end
    if isMid
    f = find(sigmask(e,:));
else
    f = find(sigmask(e,:) & ~sigmask(p,:));
    if isempty(f), f = find(sigmask(e,:)); end   % bilateral: fall back to all cluster freqs
end

    if isempty(f), continue; end
    aE(end+1,1) = mean(abs(X(e,f)));
    cE(end+1,1) = mean(abs(X(p,f)));           % identical to aE for midline, so LI = 0
end


LI = (aE - cE) ./ (aE + cE);

[~, p, ~, stats] = ttest(LI, 0);
d = stats.tstat / sqrt(stats.df + 1);
fprintf('Cluster %d: mean LI = %.3f, t(%d) = %.2f, p = %.4f, d = %.2f\n', ...
        iclust, mean(LI), stats.df, stats.tstat, p, d);