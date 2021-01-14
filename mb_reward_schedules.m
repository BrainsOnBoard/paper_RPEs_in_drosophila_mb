function R = mb_reward_schedules(ID,seed,SD,NT,NO,varargin)
%
% Generate reward schedule for td_mushroombody models
%
% Inputs
% ID - flag for the desired reward schedule
% SD - standard deviation for Gaussian white noise
% NT - # time steps
% NO - # options
%
% Outputs:
% R - reward schedule [NT x NO]

%%% Init random # stream
seeds = 10000:65000;
seeds = seeds(isprime(seeds));
rng(seeds(seed));          

if ID==0
  % Constant with noise
  amp = varargin{1};
  R = amp + SD*randn(NT,NO);
  
elseif ID==1
  % Reversal followed by equal negative followed by 2 different negatives  
  ns = floor(NT/4);
  R = zeros(NT,NO);
  R(1:ns,1) = 1;
  R(1:ns,2) = -1;
  R(ns+1:2*ns,1) = -1;
  R(ns+1:2*ns,2) = 1;
  R(2*ns+1:3*ns,1) = 0;
  R(2*ns+1:3*ns,2) = 0;
  R(3*ns+1:end,1) = 1;
  R(3*ns+1:end,2) = 0;  
  R = R + SD*randn(NT,NO);
  
elseif ID==2
  % Reversal then multiple steps (for differentiating signed-valence and signed RPE)
  ns = floor(NT/8);
  R = zeros(NT,NO);
  R(1:ns,1) = 2;
  R(1:ns,2) = -2;
  for j=2*ns:ns:8*ns
    R((j-ns)+1:j,1) = -4 + (j/ns);
    R((j-ns)+1:j,2) = 4 - (j/ns);
  end;
  if (NT-j)>0
    R(j:end,:) = repmat(R(j-1,:),NT-j,1);
  end;
  R = R + SD*randn(NT,NO);
  
elseif ID==3
  % Random walk
  amp = varargin{1};
  timescale = varargin{2};  
  R = zeros(NT,NO);    
  R = mygfilter((randn(NT,NO)),[timescale 0]);  
  R = amp * R / max(R(:)) + SD * randn(NT,NO);
  
elseif ID==4
  % Square wave
  if nargin>5, period = varargin{1};
  else period = 100; end;
  if nargin>6, offset = varargin{2};
  else offset = 0; end;
  if nargin>7, amp = varargin{3};
  else amp = 0.5; end;
  
  R = zeros(NT,NO);
  R(:,1) = offset + amp*square(2*pi*(1:NT)'/period);
  R(:,2) = offset - amp*square(2*pi*(1:NT)'/period);  
  R = R + SD*randn(NT,NO);
  
elseif ID==5
  % Blocking experiment
  if nargin>4
    amp = varargin{1};
  else
    amp = 1;
  end;
  ns = floor(NT/2);
  R = zeros(NT,NO);  
  R(1:ns,1) = amp;
  for j=2:NO
    R(:,j) = 0;
  end;
  R = R + SD*randn(NT,NO);  
  
elseif ID==6
  % Escalating reward magnitude to test linearity of reward estimate
  ns = floor(NT/2);
  R = zeros(NT,NO);
  R(1:ns,1) = 1:ns;
  R(1:ns,2) = 1:ns;
  R(ns+1:NT,1) = -1:-1:-ns;
  R(ns+1:NT,2) = -1:-1:-ns;
  R = R + SD*randn(NT,NO);     
  
elseif ID==7
  % R = [0, -1], then R = [-0.5, -1]: replicating 2013perisse.waddell paper
  % (+ve DANs needed for learning less aversive association)
  ns = floor(NT/2);
  R = zeros(NT,NO);
  R(1:ns,1) = 0;
  R(1:ns,2) = -2;
  R(ns+1:NT,1) = -1;
  R(ns+1:NT,2) = 0;
  R = R + SD*randn(NT,NO);
  
elseif ID==8
  % Sine wave
  period = 100;  
  offset = 1;
  amp = 1;
  R = zeros(NT,NO);
  R(:,1) = offset + amp*sin(2*pi*(1:NT)'/period);
  R(:,2) = offset - amp*sin(2*pi*(1:NT)'/period);  
  R = R + SD*randn(NT,NO);
  
elseif ID==9
  % Ornstein-Uhlenbeck process
  if nargin>4
    tau = varargin{1};
  else
    tau = 500;
  end;
  ns = NT;
  mu = 0;
  eta = 1/tau;  
  R = zeros(ns,NO);
  for j=2:ns 
    R(j,:) = R(j-1,:) - eta*R(j-1,:) + SD*randn(1,NO);
  end;
  R = R + mu;
  
elseif ID==10
  % Incrementing reward magnitude to demonstrate failures of old MB model
  % For Fig 2 & 3 in 2018bennett.nowotny_resc_wag
  if nargin>4
    % Exponential growth in reward (Fig 3)
    a = varargin{1};
    ns = floor(NT/2);
    R = zeros(NT,NO);
    R(1:ns,1) = exp((1:ns)*a) - 1;
    R(1:ns,2) = 0;
    R = R + SD*randn(NT,NO);
  else
    % Linear growth in reward (Fig 2)
    ns = floor(NT/2);
    R = zeros(NT,NO);
    R(1:ns,1) = (1:ns)/ns;
    R(1:ns,2) = 0;
    R = R + SD*randn(NT,NO);
  end;
  
elseif ID==11
  % Changing reward variance, constant mean  
  R = zeros(NT,NO);
  R(1:floor(NT/2),1) = 0.5;
  R(floor(NT/2)+1:end,1) = 2;
  for j=2:NO
    R(:,j) = ones(NT,1);
  end;
  R = SD * R .* randn(NT,NO);
  
elseif ID==12
  % Typical appetitive/aversive conditioning protocol (CS+, CS-, test)
  if nargin>5
    amp = varargin{1};
  else
    amp = 1;
  end;
  if mod(NT,3)~=0
    error('??? TD_MUSHROOMBODY_REWARD_SCHEDULES: ID==12, NT must be divisible by 3.')
  end;
  ns = floor(NT/3);
  R = zeros(NT,NO);  
  R(1:ns,1) = amp;
%   R(2*ns+1:3*ns,1) = amp;
  R = R + SD*randn(NT,NO);  
elseif ID==13
  % Stepped reward schedule, that moves from neutral to reward to
  % punishment and back to neutral. Reward for option 2 is much lower, 
  % encouraging only option 1 to be chosen (ideally, set decision=1 in the 
  % simulation script)
  ns = floor(NT/9);
  R = zeros(NT,NO);  
  amp = 1;
  for j=1:3
    R((j-1)*ns+1:j*ns,1) = amp*(j - 1);
  end;
  for j=4:7
    R((j-1)*ns+1:j*ns,1) = amp*(5 - j);
  end; 
  for j=8:9
    R((j-1)*ns+1:j*ns,1) = amp*(j - 9);
  end;
%   R(:,2) = -1000;
  R = R + SD*randn(NT,NO);  
elseif ID==14
  % To replicate the Real (1981) experiments (same means, different
  % variances
  amp1 = varargin{1};
  prob1 = varargin{2};
  amp2 = varargin{3};
  prob2 = varargin{4};
  R = zeros(NT,NO);  
  R(:,1) = amp1 * double(rand(NT,1) < prob1);
  R(:,2) = amp2 * double(rand(NT,1) < prob2);
elseif ID==15
  % Modification of the Real (1981) experiments (mean and SD of cue 2 are 
  % free parameters.
  amp1 = varargin{1};
  sd1 = varargin{2};
  amp2 = varargin{3};
  sd2 = varargin{4};
  R = zeros(NT,NO);
  R(:,1) = amp1 + sd1*randn(NT,1);
  R(:,2) = amp2 + sd2*randn(NT,1);
elseif ID==16
  % Reward reversals with a given probability
  stimes = varargin{1};
  amp = varargin{2};
  R = amp * ones(NT,NO);  
  R(:,2) = -amp;
  for j=2:length(stimes)
		R(stimes(j-1)+1:stimes(j),1) = amp * (-1 + 2*mod(j,2));
		R(stimes(j-1)+1:stimes(j),2) = amp * (1 - 2*mod(j,2));
	end;
	if stimes(j)<NT
		R(stimes(j)+1:end,1) = amp * (-1 + 2*mod(j+1,2));
		R(stimes(j)+1:end,2) = amp * (1 - 2*mod(j+1,2));
	end;
  R = R + SD*randn(NT,NO);
end;
