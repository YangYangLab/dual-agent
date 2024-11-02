% rebuild full sessData struct from TrialType & any other vector containing mouse data

function sdata = bm_rebuild_sessData(TrialType, otherName, otherVar)
if nargin < 3
    logging.error('rebuild_sessData requires at least 3 arguments')
    help(mfilename)

    sdata = [];
    return
end

PNS = 1; RWD = 2;

nTrials = length(TrialType);

switch lower(otherName)
case 'result'
    sdata = struct();
    sdata.n = nTrials;
    sdata.TrialType = TrialType;
    sdata.IsDeceptive = TrialType>2.5;
    sdata.Stimulus = TrialType; sdata.Stimulus(sdata.IsDeceptive) = sdata.Stimulus(sdata.IsDeceptive) - 2;
    sdata.Result = otherVar;
    sdata.Rewarded = sdata.Result==RWD;
    cm = [ 2 1 0 0; 1 2 0 0; 1 2 0 0; 2 1 0 0];
    sdata.Choice = arrayfun(@(x) cm(sdata.TrialType(x), sdata.Result(x)), 1:sdata.n);
    cm = [ 1 0; 0 1; 1 0; 0 1];
    sdata.Correctness = arrayfun(@(x) cm(sdata.TrialType(x), sdata.Choice(x)), 1:sdata.n);
    sdata.Feedback = (sdata.IsDeceptive*2) + (2-sdata.Rewarded);
    sdata.TrialSummary = (sdata.Choice==1).*(sdata.Feedback) + (sdata.Choice==2).*(4+sdata.Feedback);
case 'choice'
    rm = [RWD PNS; PNS RWD; PNS RWD; RWD PNS;];
    r = arrayfun(@(x) rm(TrialType(x), otherVar(x)), 1:nTrials);
    sdata = bm_rebuild_sessData(TrialType, 'Result', r);
case 'rewarded'
    r = other.Rewarded+1;
    sdata = bm_rebuild_sessData(TrialType, 'Result', r);
case 'correctness'
    rm = [PNS RWD; PNS RWD; RWD PNS; RWD PNS;];
    r = arrayfun(@(x) rm(TrialType(x), otherVar(x)), 1:nTrials);
    sdata = bm_rebuild_sessData(TrialType, 'Result', r);
case 'feedback'
    rm = [RWD PNS RWD PNS];
    r = (mod(otherVar, 2)==1) .* RWD + (mod(otherVar, 2)==0) .* PNS;
    sdata = bm_rebuild_sessData(TrialType, 'Result', r);
case 'trialsummary'
    rm = [RWD PNS RWD PNS RWD PNS RWD PNS];
    r = (mod(otherVar, 2)==1) .* RWD + (mod(otherVar, 2)==0) .* PNS;
    sdata = bm_rebuild_sessData(TrialType, 'Result', r);
otherwise
    logging.error('rebuild_sessData: unknown otherName')
    help(mfilename)

    sdata = [];
    return
end