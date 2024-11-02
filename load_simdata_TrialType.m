function [TrialType] = load_simdata_TrialType(pct, session)

datapath = [fileparts(which(mfilename())) './data_for_simulation'];
sessFile = fullfile(datapath, sprintf('TrialType-%2d-S%02d.mat', pct, session));
TrialType = load(sessFile);
TrialType = TrialType.trialtype;
