function out = mb_mv_b(gamma,seed,nt,rs_flag,epskm,no)
%
% Mixed valence MB model. Same as mb_mv_a, but modified so that predictions 
% are updated for all cues on each trial in order to examine contributions
% to reward predictions from each single KC. There is no decision making in
% this verion of the model.
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

%%% Init random # stream
seeds = 10000:65000;
seeds = seeds(isprime(seeds));
rng(seeds(seed));          

%%%% Reward schedule
r = rs_flag;

%%% Network setup
nk = no * 10; % # KCs
sparseness = 1 / no; % KC sparseness
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

%%% Run simulation
for j=1:nt  % Loop over trials
  % Compute MBON firing rates
  for stim=1:no    
    map(j,stim) = wkmap(:,:,j) * s(:,stim);
    mav(j,stim) = wkmav(:,:,j) * s(:,stim);        
        
    % Positive value    
    go(j,stim) = map(j,stim);
    % Negative value
    nogo(j,stim) = mav(j,stim);
  end;        
  
  % Compute DAN firing rates
  for stim=1:no
    dap(j,stim) = max(0,wkdap * s(:,stim) - wmapdap * map(j,stim) + wmavdap * mav(j,stim) + r(j,stim));
    dav(j,stim) = max(0,wkdav * s(:,stim) - wmavdav * mav(j,stim) + wmapdav * map(j,stim) - r(j,stim));
  end;
      
  % Update KC->MBON weights (except on last trial)
  if j<nt
    dwkmap = zeros(1,nk);
    dwkmav = zeros(1,nk);
    for stim=1:no      
      dwkmap = dwkmap + epskm * s(:,stim)' .* (dap(j,stim) - dav(j,stim));
      dwkmav = dwkmav + epskm * s(:,stim)' .* (dav(j,stim) - dap(j,stim));
    end;
    wkmap(:,:,j+1) = max(0,wkmap(:,:,j) + dwkmap);
    wkmav(:,:,j+1) = max(0,wkmav(:,:,j) + dwkmav);
  end;
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
out.r = r;
out.s = s;
