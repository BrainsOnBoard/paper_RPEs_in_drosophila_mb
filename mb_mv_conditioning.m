function out = mb_mv_conditioning(seed,rewards,epskm,varargin)
%
% Model as in mb_mv_a, for the specific case of simulating Drosophila
% conditioning experiments: 1) train CS+, 2) train CS-, 3) test
% intervene_id - choose a cell type to block:
%          1 - map
%          2 - mav
%          3 - dap
%          4 - dav

%%% Init random # stream
seeds = 10000:65000;
seeds = seeds(isprime(seeds));
rng(seeds(seed));          

%%% Model parameters
nk = 100;
no = 2;
sparseness = 0.1;
nt = 30;
gamma = 1;

%%% Process options
intervene_id = 0; % Default: no intervention
pr = 'eq7'; % Default: Eq. 7 learning rule
if nargin>3
  j = 1;
  while j<=numel(varargin)
    if strcmp(varargin{j},'intervene_id')
      intervene_id = varargin{j+1}; % Which cell type
      intrvn_type = varargin{j+2}; % 1: multiplicative; 0: additive
      if isscalar(varargin{j+3})
        intrvn_strength = ones(nt,1)*varargin{j+3}; % Strength of intervention
      else
        intrvn_strength = varargin{j+3}; % Strength of intervention
      end;
      j = j + 4;
    elseif strcmp(varargin{j},'plasticity_rule')
       pr = varargin{j+1}; j = j + 2; % Which learning rule to use
    end;
  end;
end;

%%%% Reward schedules
r = rewards;

% Init synaptic weights
wkmap = zeros(1,nk,nt); % KC->M+
wkmav = zeros(1,nk,nt); % KC->M-
wkmap(:,:,1) = 0.1*rand(1,nk);
wkmav(:,:,1) = 0.1*rand(1,nk);
wkdap = gamma * ones(1,nk);  % KC->D+
wkdav = gamma * ones(1,nk);  % KC->D-
wmapdap = 1; % M+->D+
wmavdap = 1; % M-->D+
wmapdav = 1; % M+->D-
wmavdav = 1; % M-->D-

% Softmax temperature
T = 0.2;
beta = 1 / T;

%%% Generate KC responses to cues
s = zeros(nk,no);
for j=1:no
  s(:,j) = double(rand(nk,1)<sparseness);
  s(:,j) = s(:,j) / sum(s(:,j)) * 10;    
end;
%%% Allocate memory
dap = zeros(nt,1);
dav = zeros(nt,1);
map = zeros(nt,no);
mav = zeros(nt,no);
go = zeros(nt,no);
nogo = zeros(nt,no);
decision = zeros(nt,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%% Train CS+ (stim 1)
%%%
for j=1:(nt/3)
  % Compute MBON firing rates
  stim = 1;
  map(j,stim) = wkmap(:,:,j) * s(:,stim);
%   mav(j,stim) = wkmav(:,:,j) * s(:,stim);
  mav(j,stim) = max(0,wkmav(:,:,j) * s(:,stim) - map(j,stim));
  
  % Apply intervention to MBONs (if any)
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
  % Reward prediction
  mdiff(stim) = go(j,stim) - nogo(j,stim);
  
  % Decision using softmax
  decision(j) = stim;
  
  % DAN firing rates
  dap(j) = max(0,wkdap * s(:,decision(j)) - wmapdap * map(j,decision(j)) + wmavdap * mav(j,decision(j)) + r(j,decision(j)));
  dav(j) = max(0,wkdav * s(:,decision(j)) - wmavdav * mav(j,decision(j)) + wmapdav * map(j,decision(j)) - r(j,decision(j)));   
  
  % Apply intervention to DANs (if any)
  if any(intervene_id==3)
    dap(j) = max(0,intrvn_type * dap(j) * intrvn_strength(j) + (1-intrvn_type) * (dap(j) + intrvn_strength(j)));
  end;
  if any(intervene_id==4)
    dav(j) = max(0,intrvn_type * dav(j) * intrvn_strength(j) + (1-intrvn_type) * (dav(j) + intrvn_strength(j)));
  end;
  
  % Update KC->MBON synaptic weights
  if strcmp(pr,'eq7')
    wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
    wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
  elseif strcmp(pr,'eq8')
    wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (wkdav * s(:,decision(j)) - dav(j)));
    wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (wkdap * s(:,decision(j)) - dap(j)));
  end;
end;

%%%
%%% Train CS- (stim 2)
%%%
for j=(nt/3 + 1):(2*nt/3)
  % Compute MBON firing rates
  stim = 2;
  map(j,stim) = wkmap(:,:,j) * s(:,stim);
  mav(j,stim) = wkmav(:,:,j) * s(:,stim);
  
  % Apply intervention to MBONs (if any)
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
  % Reward prediction
  mdiff(stim) = go(j,stim) - nogo(j,stim);
  
  % Decision using softmax
  decision(j) = stim;
  
  % DAN firing rates
  dap(j) = max(0,wkdap * s(:,decision(j)) - wmapdap * map(j,decision(j)) + wmavdap * mav(j,decision(j)) + r(j,decision(j)));
  dav(j) = max(0,wkdav * s(:,decision(j)) - wmavdav * mav(j,decision(j)) + wmapdav * map(j,decision(j)) - r(j,decision(j)));   
  
  % Apply intervention to DANs (if any)
  if any(intervene_id==3)
    dap(j) = max(0,intrvn_type * dap(j) * intrvn_strength(j) + (1-intrvn_type) * (dap(j) + intrvn_strength(j)));
  end;
  if any(intervene_id==4)
    dav(j) = max(0,intrvn_type * dav(j) * intrvn_strength(j) + (1-intrvn_type) * (dav(j) + intrvn_strength(j)));
  end;  
  
  % Update KC->MBON synaptic weights
if strcmp(pr,'eq7')
    wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
    wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
  elseif strcmp(pr,'eq8')
    wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (wkdav * s(:,decision(j)) - dav(j)));
    wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (wkdap * s(:,decision(j)) - dap(j)));
  end;
end;

%%%
%%% Test
%%%
for j=(2*nt/3+1):nt
  % Compute MBON firing rates for all cues
  for stim=1:no
    map(j,stim) = wkmap(:,:,j) * s(:,stim);
    mav(j,stim) = wkmav(:,:,j) * s(:,stim);        
    
    % Apply intervention to MBONs (if any)
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
    % Reward prediction
    mdiff(stim) = go(j,stim) - nogo(j,stim);        
  end;        
  
  % Decision using softmax
  decision(j) = double(rand < exp(mdiff(2)*beta)/(exp(mdiff(1)*beta)+exp(mdiff(2)*beta))) + 1;  
  
  % DAN firing rates
  dap(j) = max(0,wkdap * s(:,decision(j)) - wmapdap * map(j,decision(j)) + wmavdap * mav(j,decision(j)) + r(j,decision(j)));
  dav(j) = max(0,wkdav * s(:,decision(j)) - wmavdav * mav(j,decision(j)) + wmapdav * map(j,decision(j)) - r(j,decision(j)));   
  
  % Apply intervention to DANs (if any)
  if any(intervene_id==3)
    dap(j) = max(0,intrvn_type * dap(j) * intrvn_strength(j) + (1-intrvn_type) * (dap(j) + intrvn_strength(j)));
  end;
  if any(intervene_id==4)
    dav(j) = max(0,intrvn_type * dav(j) * intrvn_strength(j) + (1-intrvn_type) * (dav(j) + intrvn_strength(j)));
  end;
  
  % Update KC->MBON synaptic weights (except on last trial)
  if j<nt
    if strcmp(pr,'eq7')
      wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
      wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
    elseif strcmp(pr,'eq8')
      wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (wkdav * s(:,decision(j)) - dav(j)));
      wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (wkdap * s(:,decision(j)) - dap(j)));
    end;
  end;
end;

% Create output struct
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
