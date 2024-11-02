% perform multiple runs of the cognitive model

clear all

iteration = 50;
pct = 80; nTrials = 8000; nTrialsInSession = 500; Optogenetics = [0]; savename = '80%';
% pct = 90; nTrials = 6000; nTrialsInSession = 500; Optogenetics = [0]; 
% pct = 95; nTrials = 4000; nTrialsInSession = 500; Optogenetics = [0]; 
% pct = 80; nTrials = 8000; nTrialsInSession = 500; Optogenetics = [1 5]; savename = 'optoStim';

%read 8000 trial sequence from real data
TrialType = [];
for k = 1:nTrials/nTrialsInSession
    trials = load_simdata_TrialType(pct, k);
    TrialType = [TrialType trials(1:nTrialsInSession)];
end

rng(1); % for reproducibility

models = cell(1, iteration);
%alpha = [0.019 0.028 0.003];
for iter = 1:iteration
    onerun = bm_cognitive_ProbLearn_multiagent(TrialType, 'optogenetics', Optogenetics, 'noise', 0.04); %, 'alpha', alpha);
    models{iter} = onerun;
end
