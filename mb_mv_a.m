function out = mb_mv_a(gamma,seed,nt,rs_flag,epskm,varargin)
%
% Mixed valence MB model.
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
intervene_id = 0;
choose1 = false;
no = 2;
period = 100;
nk = 20; % # KCs
lambda = 0;
pr = 'eq7';
%%% Update default parameters with custom options
if nargin>5
  j = 1;
  while j<=numel(varargin)
    if strcmp(varargin{j},'intervene_id')
      % For simulating genetic interventions (e.g. shibire/dTrpA1)
      intervene_id = varargin{j+1}; % Which cell type
      intrvn_type = varargin{j+2}; % 1: multiplicative; 0: additive
      if isscalar(varargin{j+3})
        intrvn_strength = ones(nt,1)*varargin{j+3}; % Strength of intervention
      else
        intrvn_strength = varargin{j+3}; % Strength of intervention
      end;
      j = j + 4;
    elseif strcmp(varargin{j},'choose1')
      % If choose1==true, cue 1 is chosen on every trial.
      choose1 = varargin{j+1}; j = j + 2;
    elseif strcmp(varargin{j},'no')
      % Specify the number of cues
      no = varargin{j+1};
      j = j+2;
    elseif strcmp(varargin{j},'nk')
      % Specify the number of KCs
      nk = varargin{j+1};
      j = j+2;
    elseif strcmp(varargin{j},'period')
      % If using periodic reward schedule: specifiy the period
      period = varargin{j+1};
      j = j+2;
    elseif strcmp(varargin{j},'lambda')
      % If using periodic reward schedule: specifiy the period
      lambda = varargin{j+1};
      j = j+2;
    elseif strcmp(varargin{j},'plasticity_rule')
       pr = varargin{j+1}; j = j + 2; % Which learning rule to use
    else
      error('???MB_MV_A: Optional arguments not recognised.');
    end;
  end;
end;

%%% Init random # stream
seeds = 10000:65000;
seeds = seeds(isprime(seeds));
rng(seeds(seed));          

%%%% Reward schedules..............(ID,SD,NT,NO)
if isscalar(rs_flag) % if generating a new reward schedule
  if rs_flag==4
    r = mb_reward_schedules(rs_flag,seed,0.1,nt,2,period);
  else
    r = mb_reward_schedules(rs_flag,seed,0.1,nt,2);
  end;
else % if using a pregenerated reward schedule
  r = rs_flag;
end;

%%% Network setup
sparseness = 0.5; % KC sparseness
% Softmax temperature
T = 0.2;
beta = 1 / T;

%%% Initialise synaptic weights
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

%%% Generate KC responses to cues
s = zeros(nk,no);
for j=1:no
%   s(:,j) = double(rand(nk,1)<sparseness);
  s(floor((j-1)*sparseness*nk)+1:floor(j*sparseness*nk),j) = 1;
  s(:,j) = s(:,j) / sum(s(:,j)) * 10;    
end;

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
for j=1:nt  % Loop over trials
  % Compute MBON firing rates and reward predictions (mdiff)
  for stim=1:no
    map(j,stim) = wkmap(:,:,j) * s(:,stim);
    mav(j,stim) = wkmav(:,:,j) * s(:,stim);
        
    % For "genetic" interventions
    if any(intervene_id==1)
      map(j,stim) = max(0,intrvn_type * map(j,stim) * intrvn_strength(j) + (1-intrvn_type) * (map(j,stim) + intrvn_strength(j)));
    end;
    if any(intervene_id==2)
      mav(j,stim) = max(0,intrvn_type * mav(j,stim) * intrvn_strength(j) + (1-intrvn_type) * (mav(j,stim) + intrvn_strength(j)));
    end;
    
    % Positive value    
    go(j,stim) = map(j,stim);
    % Negative value
    nogo(j,stim) = mav(j,stim);      
  end;        
  
  % Make decision
  if choose1
    decision(j) = 1;
  else
    for stim=1:no
      probs(stim) = exp((go(j,stim)-nogo(j,stim))*beta) / sum(exp((go(j,:)-nogo(j,:))*beta));
    end;
    probs = probs / sum(probs); % Ensure probs is normalised to 1 (to avoid rounding errors)
    randchoice = rand;
    flag = 1; k = 1;
    while flag
      if k>no % FOR DEBUGGING
        error('???MB_MV_A: no choice made.');
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
  
  % For "genetic" interventions
  if any(intervene_id==3)
    dap(j) = max(0,intrvn_type * dap(j) * intrvn_strength(j) + (1-intrvn_type) * (dap(j) + intrvn_strength(j)));
  end;
  if any(intervene_id==4)
    dav(j) = max(0,intrvn_type * dav(j) * intrvn_strength(j) + (1-intrvn_type) * (dav(j) + intrvn_strength(j)));
  end;
  
  % Update KC->MBON weights (except on last trial)
  if j<nt
    if strcmp(pr,'eq7')
      wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
      wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
    elseif strcmp(pr,'eq8')
      wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (wkdav * s(:,decision(j)) - dav(j)));
      wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (wkdap * s(:,decision(j)) - dap(j)));
    end;
  end;
  
  % Update summed reward and RPE
  sr = sr + r(j,decision(j));
end;

%%% Creater output struct
out.map = map;
out.mav = mav;
out.dap = dap;
out.dav = dav;
out.go = go;
out.nogo = nogo;
out.wkmap = wkmap;
out.wkmav = wkmav;
out.decision = decision;
n1 = decision==1; n2 = decision==2; 
out.pi = (n1-n2) ./ (n1+n2);
out.r = r;
out.s = s;
out.sr = sr;
