function out = mb_mv_blocking(seed,epskm,varargin)
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
nk = 40;
no = 2;
sparseness = 0.5;
nt = 30;
gamma = 1;
pflip = [0 0];

%%% Process options
intervene_id = 0;
if nargin>2
  elseif strcmp(varargin{1},'pflip')
    pflip = varargin{2};    
  elseif strcmp(varargin{1},'intervene_id')
    intervene_id = varargin{2}; % Which cell type
    intrvn_type = varargin{3}; % 1: multiplicative; 0: additive
    if isscalar(varargin{4})
      intrvn_strength = ones(nt,1)*varargin{4}; % Strength of intervention
    else
      intrvn_strength = varargin{4}; % Strength of intervention
    end;
  end;    
end;

%%%% Reward schedules
r = zeros(nt,no);
r(1:10,1) = 1;
r(11:20,2) = 1;
r = r + 0.1*randn(30,2);

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
ss = zeros(nk,no);
for j=1:no
  ss((j-1)*nk/2+1:j*nk/2,j) = double(rand(nk/2,1)<sparseness);
%   s(floor((j-1)*sparseness*nk)+1:floor(j*sparseness*nk),j) = 1;
  s(:,j) = ss(:,j) / sum(ss(:,j)) * 10;  
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
  wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
  wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
end;

%%%
%%% Train CS- (stim 2)
%%%
flipind1 = [rand(nk/2,1) < pflip(1); rand(nk/2,1) < pflip(2)];
sss = sum(ss,2); %sss(~sparseind) = 0; 
flipind2 = sss>0 & flipind1;
flipind3 = find(~(sss>0) & flipind1);
sss(flipind2) = 1 - sss(flipind2);
if any(flipind2)
  ind = randperm(length(flipind3));
  flipind4 = flipind3(ind(1:min(length(flipind3),sum(double(flipind2)))));  
  sss(flipind4) = 1 - sss(flipind4);
end;
s(:,2) = sss / sum(sss) * 20;

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
  wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
  wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
end;

%%%
%%% Test
%%%
for j=1:no
  s(:,j) = ss(:,j) / sum(ss(:,j)) * 10;
end;
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
  decision(j) = double(rand < exp(mdiff(2)*beta)/(1+exp(mdiff(2)*beta))) + 1;  
  
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
    wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + epskm * s(:,decision(j))' .* (dap(j) - dav(j)));
    wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + epskm * s(:,decision(j))' .* (dav(j) - dap(j)));
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
out.sss = sss/sum(sss)*20;
