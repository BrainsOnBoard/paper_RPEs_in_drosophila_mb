function out = mb_mv_c(gamma,seed,nt,rs_flag,varargin)
%
% Same as mb_mv_a, but modified so that predictions for a 
% given cue are updated to be equal to the reward just received.
%
% Inputs:
%    gamma - KC->DAN synaptic weight
%     seed - an integer, N, that selects a prime number to seed the random
%            number generator
%       nt - # trials
%  rs_flag - reward schedule ID
%    epskm - learning rate
%
% Outputs:
%  out - struct containing numerous fields (see bottom of script)

%%% Set defaults for optional parameters
choose1 = false;
no = 2;

%%% Update default parameters with custom options
if nargin>4
  j = 1;
  while j<=numel(varargin)    
    if strcmp(varargin{j},'choose1')
      % If choose1==true, cue 1 is chosen on every trial.
      choose1 = varargin{j+1}; j = j + 2;
    elseif strcmp(varargin{j},'no')
      % Specify the number of cues
      no = varargin{j+1};
      j = j+2;
    elseif strcmp(varargin{j},'period')
      % If using periodic reward schedule: specifiy the period
      period = varargin{j+1};
      j = j+2;
    else
      error('???MB_MV_C: Optional arguments not recognised.');
    end;
  end;
end;

%%% Init random # stream
seeds = 10000:65000;
seeds = seeds(isprime(seeds));
rng(seeds(seed));          

%%%% Reward schedule
r = rs_flag;

%%% Network setup
nk = no * 100;
sparseness = 1 / no;

% Initialise synaptic weights
wkmap = zeros(1,nk,nt); % KC -> M+
wkmav = zeros(1,nk,nt); % KC -> M-
wkmap(:,:,1) = 0.1*rand(1,nk);
wkmav(:,:,1) = 0.1*rand(1,nk);
wkdap = gamma * ones(1,nk); % KC -> D+
wkdav = gamma * ones(1,nk); % KC -> D-
wmapdap = 1; % M+ -> D+
wmavdap = 1; % M- -> D+
wmapdav = 1; % M+ -> D-
wmavdav = 1; % M- -> D-
% Softmax temperature
T = 0.2;
beta = 1 / T;

%%% Generate KC responses to cues
s = zeros(nk,no);
for j=1:no
  s(floor((j-1)*sparseness*nk)+1:floor(j*sparseness*nk),j) = 1;
  s(:,j) = s(:,j) / sum(s(:,j)) * 10;    
end;
sums = sum(s,1);

%%% Allocate memory for firing rates
dap = zeros(nt,1);
dav = zeros(nt,1);
map = zeros(nt,no);
mav = zeros(nt,no);
go = zeros(nt,no);
nogo = zeros(nt,no);
decision = zeros(nt,1);
sr = 0;
probs = zeros(no,1);

%%% Run simulation
for j=1:nt  
  % Compute MBON firing rates and reward predictions (mdiff) for all options
  for stim=1:no
    map(j,stim) = wkmap(:,:,j) * s(:,stim);
    mav(j,stim) = wkmav(:,:,j) * s(:,stim);        
            
    % Positive value    
    go(j,stim) = map(j,stim);
    % Negative value
    nogo(j,stim) = mav(j,stim);
  end;        
  
  % Decision using softmax  
  if choose1
    decision(j) = 1;
  else
    for stim=1:no
      probs(stim) = exp((go(j,stim)-nogo(j,stim))*beta) / sum(exp((go(j,:)-nogo(j,:)).*beta));
    end;
    probs = probs / sum(probs); % Ensure probs is normalised to 1 (to avoid rounding errors)
    randchoice = rand;
    flag = 1; k = 1;
    while flag
      if k>no % FOR DEBUGGING
        error('???MB_MV_C: no choice made.');
%         keyboard;
%         fprintf('K ABOVE NO\n');
%         fprintf('sumprobs = %.5f     randchoice = %f\n',sum(probs),randchoice);
%         fprintf('no = %d      seed = %d\n',no,seed);
%         for zz=1:nk
%           fprintf('%f\n',s(zz,17));
%         end;
      end;
      if randchoice<sum(probs(1:k))
        decision(j) = k;
        flag = 0;
      end;
      k = k + 1;
    end;
  end;
  
  % Compute DAN firing rates
  dap(j) = max(0,wkdap * s(:,decision(j)) - wmapdap * map(j,decision(j)) + wmavdap * mav(j,decision(j)) + r(j,decision(j)));
  dav(j) = max(0,wkdav * s(:,decision(j)) - wmavdav * mav(j,decision(j)) + wmapdav * map(j,decision(j)) - r(j,decision(j)));   

  % Update KC->MBON weights (except on last trial)
  if j<nt
    ind = s(:,decision(j)) > 0;
    wkmap(:,:,j+1) = wkmap(:,:,j);
    wkmav(:,:,j+1) = wkmav(:,:,j);
    wkmap(:,ind,j+1) = max(0,r(j,decision(j))) / sums(decision(j));
    wkmav(:,ind,j+1) = -min(0,r(j,decision(j))) / sums(decision(j));    
  end;
  
  % Update summed reward and RPE
  sr = sr + r(j,decision(j));
end;

%%% Create output struct
out.map = map;
out.mav = mav;
out.dap = dap;
out.dav = dav;
out.go = go;
out.nogo = nogo;
out.wkmap = wkmap;
out.wkmav = wkmav;
out.decision = decision;
out.r = r;
out.s = s;
out.sr = sr;
