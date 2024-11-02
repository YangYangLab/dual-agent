% a cognitive model, combining rsWSLS and model-based RL

function [model] = bm_cognitive_ProbLearn_v2(TrialType, rsWSLS, RL, params, opts)
arguments
    TrialType           = [];

    rsWSLS.ratePE    = [0.190 0.135];
    rsWSLS.surprise  = 0.65;     % surprise threshold
    rsWSLS.initWiSt  = 0.65;
    rsWSLS.initLoSh  = 1-0.45;
    rsWSLS.thetaWin  = 0.050;
    rsWSLS.thetaLose = 0.035;
    
    RL.alpha            = [0.014 0.031 0.003];
    RL.noise            = 0.04;
    RL.spontaneous      = 1e-3;
    RL.temperature      = 0.12;
    RL.diffusion        = 0.00005;
    
    params.perseverance = 0.60;
    params.exploration  = 0.05;
    params.optogenetics = 0;%[ 1 2 4 ];
    opts.Plot           = false;
    opts.Debug          = false;

end

%% Generate trials if not provided

if isempty(TrialType)
    nTrials = 8000; nStimulus = 2;
    Stimulus = randi(nStimulus, [1 nTrials]); 
    IsDeceptive = rand(1,nTrials) >= 0.895; % rarity level
    TrialType = Stimulus; TrialType(IsDeceptive) = TrialType(IsDeceptive) + nStimulus;
    Optogenetics = IsDeceptive;
else
    nTrials = length(TrialType); nStimulus = max(TrialType);
    IsDeceptive = TrialType > (nStimulus/2+0.1);
    Stimulus = TrialType; Stimulus(IsDeceptive) = Stimulus(IsDeceptive) - nStimulus/2;
    Optogenetics = IsDeceptive;
end

%% Initialize model values

% --- shared definitions ---

    LEFT = 1; RIGHT = 2;
    REWARD = 1; PUNISH = 0;
    GetCommonFeedback = [REWARD PUNISH; PUNISH REWARD; REWARD PUNISH; PUNISH REWARD; ];
    GetRareFeedback = 1-GetCommonFeedback;

% --- System 1: rsWSLS model ---

    STAY = 1; SWITCH = -1;

    value = rand(nStimulus, 2)*1e-3;
    saliency = ones(1, nStimulus); 
    gamma = rsWSLS.ratePE(1);

    pWiSt = rsWSLS.initWiSt;
    pLoSh = rsWSLS.initLoSh;
    thetaW = rsWSLS.thetaWin;
    thetaL = rsWSLS.thetaLose;
    finalWiSt = rsWSLS.initWiSt;
    finalLoSh = rsWSLS.initLoSh;
    deltaWiSt = 0; deltaLoSh = 0; 
    surprisePE = rsWSLS.surprise;

    memChoice   = nan(1, nStimulus); 
    memOutcome  = nan(1, nStimulus);
    memRare     = nan(1, nStimulus);

% --- System 2: Q-Learning ---

    potentiation = RL.alpha(1);
    depressionRw = RL.alpha(2);
    depressionPn = RL.alpha(3);
    
    readoutNoise = RL.noise;
    spontaneous  = RL.spontaneous;
    temperature  = RL.temperature;
    diffuse      = RL.diffusion;

    Q = abs(randn(2,2).*1e-4);

% --- general & meta-cognition system ---

    saturation = 1 - params.exploration;
    perseverance = params.perseverance;

%% tracker

    tracker = InitializeTracker(...
        'Optogenetics', 'Choice', 'Outcome', 'InternalRarityState', ...
        'value', 'pWiSt', 'pLoSh', 'deltaWiSt', 'deltaLoSh', 'finalWiSt', 'finalLoSh', ...
        'Q', ...
        'dw', 'Choice1', 'Choice2', 'Choice3');

%% perform choices

for t = 1:nTrials

    % --- get stimulus and optogenetic info

    s = Stimulus(t);
    opto = Optogenetics(t) * params.optogenetics;

    % --- memory recall 

    LastChoice  = memChoice(s);
    LastOutcome = memOutcome(s);
    
    % --- rs-WSLS
    if any(opto==1), surprisePE = nan; else, surprisePE = rsWSLS.surprise; end
    if any(opto==2), LastRarityState = nan; else, LastRarityState = memRare(s);  end
    
    dice = rand();
    if LastRarityState == 0

        if LastOutcome == REWARD 
            Decision = (dice<=pWiSt) * STAY + (dice>pWiSt) * SWITCH;
        elseif LastOutcome == PUNISH 
            Decision = (dice<=pLoSh) * SWITCH + (dice>pLoSh) * STAY;
        else 
            Decision = nan; 
        end

    elseif LastRarityState == 1

        pWiStRare = pWiSt + deltaWiSt;
        pLoShRare = pLoSh + deltaLoSh;

        if LastOutcome == REWARD 
            Decision = (dice<=pWiStRare) * STAY + (dice>pWiStRare) * SWITCH;
        elseif LastOutcome == PUNISH 
            Decision = (dice<=pLoShRare) * SWITCH + (dice>pLoShRare) * STAY;
        else 
            Decision = nan; 
        end
    
    elseif isnan(LastRarityState)

        Decision = (dice<=perseverance) * STAY + (dice>perseverance) * SWITCH; 

    end

    switch Decision
        case STAY   , Choice1 = LastChoice;
        case SWITCH , Choice1 = 3-LastChoice;
        otherwise,    Choice1 = nan;
    end
    if isnan(Choice1), Choice1 = randi(2, [1 1]); end
    
    if any(opto==3), Choice1 = randi(2, [1 1]); end
    
    % --- Q-Learning
    
    if any(opto==4) 

        Choice2 = randi(2); 
        dwRL = 0;
        NoRL = true;

    else

        normQ = exp(Q(s,:)/temperature) / sum(exp(Q(s,:)/temperature));
        activation = normQ(1) - normQ(2) + normrnd(0, readoutNoise);
        if activation>0,   Choice2 = LEFT;
        else,              Choice2 = RIGHT; end
    
        entropy = -sum(normQ .* log2(normQ));
        dwRL = entropy;

        NoRL = false;

    end
    

    % --- MAKE A FINAL CHOICE

    Choice3 = randi(2);
    exploration = params.exploration;
    if any(opto==7), exploration = 1; end

    dw = [dwRL*(1-exploration) (1-dwRL)*(1-exploration) exploration];
    Choice = randsample([Choice1, Choice2, Choice3], 1, true, dw);

    % --- outcome

    if ~IsDeceptive(t)
        Outcome = GetCommonFeedback(s, Choice);
    else
        Outcome = GetRareFeedback(s, Choice);
    end

    % --- rsWSLS 
    FastExpect = value(s, Choice) * saliency(s);
    PE = Outcome - FastExpect;
    value(s, Choice) = value(s, Choice) + gamma * PE * saliency(s);

    if isnan(surprisePE)
        InternalRarityState = nan; 
    else
        InternalRarityState = double(abs(PE) >= surprisePE);
    end

    if InternalRarityState == 0
        finalWiSt = min(0.65 + PE, saturation);
        finalLoSh = 1 - (perseverance + 0.02 * log1p(100 * abs(PE)));
    elseif InternalRarityState == 1
        deltaWiSt = 0.65 * abs(PE) - 0.39;
        deltaLoSh = 0.02 * abs(PE)^(-0.7);
    end
    
    if Outcome == REWARD
        pWiSt = pWiSt + thetaW * (finalWiSt - pWiSt);        
    else
        pLoSh = pLoSh + thetaL * (finalLoSh - pLoSh);
    end

    % --- Q-Learning
    if ~NoRL
        R = Outcome;
        SC = randg() * spontaneous;
        
        if     Outcome == REWARD && Choice==LEFT 
            Q(s,LEFT)  = Q(s,LEFT)  + potentiation .* Q(s,LEFT) .* (R  - Q(s,LEFT) );
            Q(s,RIGHT) = Q(s,RIGHT) + depressionRw .*         1 .* (SC - Q(s,RIGHT));
        elseif Outcome == REWARD && Choice==RIGHT
            Q(s,RIGHT) = Q(s,RIGHT) + potentiation .* Q(s,RIGHT).* (R  - Q(s, RIGHT));        
            Q(s,LEFT)  = Q(s,LEFT)  + depressionRw .*         1 .* (SC - Q(s, LEFT));
        else
            Q(s,LEFT)  = Q(s,LEFT)  + depressionPn .*         1 .* (SC - Q(s, LEFT));
            Q(s,RIGHT) = Q(s,RIGHT) + depressionPn .*         1 .* (SC - Q(s, RIGHT));
        end
    
    end

    % --- update memory
    memChoice(s)  = Choice;
    memOutcome(s) = Outcome;
    memRare(s)    = InternalRarityState;

    % --- track this trial

    tracker = AddTrialData(tracker, t, ...
        'Optogenetics', opto, ...
        'Choice', Choice, 'Outcome', Outcome, 'InternalRarityState', InternalRarityState, ...
        'value', value, ...
        'pWiSt', pWiSt, 'pLoSh', pLoSh, 'deltaWiSt', deltaWiSt, 'deltaLoSh', deltaLoSh, ...
        'finalWiSt', finalWiSt, 'finalLoSh', finalLoSh, ...
        'Q', Q, 'uncertainty', dwRL, 'dw', dw, ...
        'Choice1', Choice1, 'Choice2', Choice2, 'Choice3', Choice3 ...
        );

end    

model.run = tracker;
model.beh = bm_rebuild_sessData(TrialType, 'Choice', tracker.Choice);

end

function tracker = InitializeTracker(varargin)
    tracker = struct();
    for k = 1:nargin
        tracker.(varargin{k}) = [];
    end
end

function tracker = AddTrialData(tracker, t, varargin)

    for k = 1:2:nargin-2
        name = varargin{k};
        value= varargin{k+1};

        if      isscalar(value)
            tracker.(name)(t) = value;
        elseif  isvector(value)
            if ~isrow(value), value = value'; end
            tracker.(name)(t,:) = value;
        elseif ismatrix(value)
            tracker.(name) = cat(3, tracker.(name), value);
        end

    end
    
end
