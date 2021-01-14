function out = mb_vs(gamma,seed,nt,rs_flag,epskm,varargin)
%
% Valence specific MB model.
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
flag_reward_sign = 1;
flag_plasticity_rule = 1;
intervene_id = 0;
choose1 = false;
lambda = 11.5;
no = 2;
nk = 10*no;

%%% Update default parameters with custom options
if nargin>5
  j = 1;
  while j<=numel(varargin)
    if strcmp(varargin{j},'rewardsign')
      % 1 - excitatory reward signals: e.g. R+ excites D+, R- excites D-
      % 0 - inhibitory reward signals: e.g. R+ % inhibits D-, R- inhibits D+
      flag_reward_sign = varargin{j+1}; j = j + 2;    
    elseif strcmp(varargin{j},'plasticity_rule')
      % 1 - delta_w+ ~ (lambda - d-)
      % 0 - delta_w+ ~ (w_k*k - d-)
      flag_plasticity_rule = varargin{j+1}; j = j + 2;    
    elseif strcmp(varargin{j},'intervene_id')
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
    elseif strcmp(varargin{j},'lambda')
      % Set the constant term in the modified plasticity rule, P^VSlambda
      lambda = varargin{j+1}; j = j + 2;
    elseif strcmp(varargin{j},'no')
      % Set the number of cues
      no = varargin{j+1}; j = j + 2;
    elseif strcmp(varargin{j},'nk')
      % Set the number of cues
      nk = varargin{j+1}; j = j + 2;
    else
      error('???MB_VS: Optional arguments not recognised.');
    end;    
  end;
end;

%%% Seed random number generator
seeds = 10000:65000;
seeds = seeds(isprime(seeds));
rng(seeds(seed));

%%% Model parameters
sparseness = 1 / no;0.5; % KC sparseness
% Softmax temperature
T = 0.2;
beta = 1 / T;

%%% Reward schedules
if isscalar(rs_flag) % if generating a new reward schedule
  r = mb_reward_schedules(rs_flag,seed,0.1,nt,2);
else % if using a pregenerated reward schedule
  r = rs_flag;
end;

%%% Initialise synaptic weights
wkmap = zeros(1,nk,nt); % KC -> M+
wkmav = zeros(1,nk,nt); % KC -> M-
wkmap(:,:,1) = 0.1*rand(1,nk); 
wkmav(:,:,1) = 0.1*rand(1,nk);
wkdap = gamma*ones(1,nk); % KC -> D+
wkdav = gamma*ones(1,nk); % KC -> D-
wmdap = 1; % M+ -> D-
wmdav = 1; % M- -> D+

%%% Generate KC responses to cues
s = zeros(nk,no);
for j=1:no
%   s(:,j) = double(rand(nk,1)<sparseness);
  s(floor((j-1)*sparseness*nk)+1:floor(j*sparseness*nk),j) = 1;
%   s(:,j) = s(:,j) / sum(s(:,j)) * sparseness * nk;  
  s(:,j) = s(:,j) / sum(s(:,j)) * 10;  
end;

%%% Allocate memory for firing rates
dap = zeros(nt-1,1);
dav = zeros(nt-1,1);
map = zeros(nt-1,2);
mav = zeros(nt-1,2);
go = zeros(nt-1,2);
nogo = zeros(nt-1,2);
decision = zeros(nt-1,1);
mdiff = zeros(no,1);
probs = zeros(no,1);
sr = 0;
rt = zeros(nt,1);

%%% Run simulation
for j=1:nt % Loop over trials
  % Compute MBON firing rates and reward prediction (mdiff)
  for k=1:no % Loop over cues    
    map(j,k) = wkmap(:,:,j) * s(:,k);
    mav(j,k) = wkmav(:,:,j) * s(:,k);
    
    % For "genetic" interventions
    if any(intervene_id==1)
      map(j,k) = max(0,intrvn_type * map(j,k) * intrvn_strength(j) + (1-intrvn_type) * (map(j,k) + intrvn_strength(j)));
    end;
    if any(intervene_id==2)
      mav(j,k) = max(0,intrvn_type * mav(j,k) * intrvn_strength(j) + (1-intrvn_type) * (mav(j,k) + intrvn_strength(j)));
    end;
    
    go(j,k) = map(j,k);
    nogo(j,k) = mav(j,k);
    mdiff(k) = map(j,k) - mav(j,k);    
%     if abs(mdiff(k))>2
% %       keyboard;
%       error('MB_VS: reward predictgion too large');
%     end;
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
        error('???MB_VS: no choice made.');
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
  if flag_reward_sign>0
    % For excitatory reward signals
    dap(j) = max(0,wkdap * s(:,decision(j)) + wmdap * mav(j,decision(j)) + max(0,r(j,decision(j))));
    dav(j) = max(0,wkdav * s(:,decision(j)) + wmdav * map(j,decision(j)) - min(0,r(j,decision(j))));
  else
    % For inhibitory reward signals
    dap(j) = max(0,wkdap * s(:,decision(j)) + wmdap * mav(j,decision(j)) + min(0,r(j,decision(j))));
    dav(j) = max(0,wkdav * s(:,decision(j)) + wmdav * map(j,decision(j)) - max(0,r(j,decision(j))));
  end;
  
  % For "genetic" interventions
  if any(intervene_id==3)
    dap(j) = max(0,intrvn_type * dap(j) * intrvn_strength(j) + (1-intrvn_type) * (dap(j) + intrvn_strength(j)));
  end;
  if any(intervene_id==4)
    dav(j) = max(0,intrvn_type * dav(j) * intrvn_strength(j) + (1-intrvn_type) * (dav(j) + intrvn_strength(j)));
  end;

  % Update KC->MBON weights (except on last trial)
  if j<nt
    if flag_plasticity_rule>0
      wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' * (lambda - dav(j)));
      wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' * (lambda - dap(j)));
    else            
      wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (wkdav*s(:,decision(j)) - dav(j)));
      wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (wkdap*s(:,decision(j)) - dap(j)));
    end;
  end;
  
  % Update summed reward and RPE
  sr = sr + r(j,decision(j));
  rt(j) = r(j,decision(j));
end;

%%% Create output struct
out.map = map;
out.mav = mav;
out.go = map;
out.nogo = mav;
out.dap = dap;
out.dav = dav;
out.wkmap = wkmap;
out.wkmav = wkmav;
out.decision = decision;
out.r = r;
out.s = s;
out.sr = sr;
out.rt = rt;