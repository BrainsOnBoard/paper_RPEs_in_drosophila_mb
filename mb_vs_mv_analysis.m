function mb_vs_mv_analysis(flag,varargin)
%
% MB_VS_MV_ANALYSIS(FLAG)
%
% Analyse data for Bennett et al. 2019.
% Valence specific (VS) and mixed valence (MV) MB models are analysed.
%
% 
%%%
%%% Output directory, where plots will be saved
%%%
fpath = '<directory name>';

if nargin>1
  % If oldflag = 1, use R2012 version of parallel computing
  % If oldflag = 0, use R2018 version of parallel computing  
  oldflag = varargin{1};
end;

%%%%
%%%% Colours for plotting
%%%%
colr1 = [.4 .7 1];
colr2 = [255 213 1]/255;
colm1 = [0 125 255]/255;
colm2 = [255 125 0]/255;
colv1 = [50 85 185]/255;
colv2 = [160 110 40]/255;
coldap = [0 185 0]/255;
coldav = [180 50 220]/255;
coldandiff = [0.3 0.3 0.3];

if flag==1
  %%%%
  %%%% Supplementary Fig. 1
  %%%% VS model cannot learn
  %%%% Learning rule: dw_x ~ k(w_k*k - d_y)
  %%%%
  
  %%% Create 10 runs
  % Reward schedule
  rew = mb_reward_schedules(13,1,0,180,2);
  % Allocate memory
  m = zeros(180,10);
  go = zeros(180,10);
  nogo = zeros(180,10);
  dap = zeros(180,10);
  dav = zeros(180,10);
  % Generate data
  for j=1:10
    q = mb_vs(1,j,45,13,2.5e-2,'rewardsign',1,'plasticity_rule',0,'choose1',true);
    go(:,j) = q.go(:,1); nogo(:,j) = q.nogo(:,1);
    dap(:,j) = q.dap; dav(:,j) = q.dav;
    m(:,j) = q.go(:,1) - q.nogo(:,1);
  end;
  
  % Plot rewards and reward predictions
  myfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  h = plot(m,'linewidth',0.2); set(h,'color',colv1);
  set(gca,'ylim',[-2.4 2.4],'xlim',[0 180],'ytick',-2:2:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  myfig(2,7.5);
  h = plot(go,'linewidth',0.2); set(h,'color',colm1);
  hold on;
  h = plot(nogo,'linewidth',0.2); set(h,'color',colm2);
  set(gca,'ylim',[0 0.5],'xlim',[0 180],'ytick',[0 0.5],'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_cue1.eps']);
  close(gcf);
  
  % Plot DAN firing rates
  myfig(2,7.5);
  h = plot(dap,'linewidth',0.2); set(h,'color',coldap);
  hold on;
  h = plot(dav,'linewidth',0.2); set(h,'color',coldav);
  set(gca,'ylim',[10 12.3],'xlim',[0 180],'ytick',10:12,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rates_cue1.eps']);
  
  % Plot DAN firing rate difference: RPE
  h = plot(dap-dav,'linewidth',0.2); set(h,'color',coldandiff);
  set(gca,'ylim',[-2.4 2.4],'xlim',[0 180],'ytick',-2:2:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rate_difference_cue1.eps']);
  close(gcf);
end;

if flag==2
  %%%%
  %%%% Fig. 2a-d
  %%%% VS model limited learning
  %%%% Learning rule: dw_x ~ k(lambda - d_y)
  %%%%
  
  %%% Create 10 runs
  % Reward schedule
  rew = mb_reward_schedules(13,1,0,180,2);
  % Allocate memory
  m = zeros(180,10);
  go = zeros(180,10);
  nogo = zeros(180,10);
  dap = zeros(180,10);
  dav = zeros(180,10);
  % Generate data
  for j=1:10
    q = mb_vs(1,j,45,13,2.5e-2,'rewardsign',1,'plasticity_rule',1,'choose1',true,'lambda',11.5);
    go(:,j) = q.go(:,1); nogo(:,j) = q.nogo(:,1);
    dap(:,j) = q.dap; dav(:,j) = q.dav;
    m(:,j) = q.go(:,1) - q.nogo(:,1);
  end;
  
  % Plot rewards and reward predictions
  myfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  h = plot([0 180],[1.5 1.5],[0 180],[-1.5 -1.5],'linewidth',0.5,'color',[0.5 0.5 0.5]);
  h = plot(m,'linewidth',0.2); set(h,'color',colv1);
  set(gca,'ylim',[-2.4 2.4],'xlim',[0 180],'ytick',-2:2:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  myfig(2,7.5);
  h = plot(go,'linewidth',0.2); set(h,'color',colm1);
  hold on;
  h = plot(nogo,'linewidth',0.2); set(h,'color',colm2);
  set(gca,'ylim',[0 2],'xlim',[0 180],'ytick',0:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_cue1.eps']);
  close(gcf);
  
  % Plot DAN firing rates
  myfig(2,7.5);
  h = plot(dap,'linewidth',0.2); set(h,'color',coldap);
  hold on;
  h = plot(dav,'linewidth',0.2); set(h,'color',coldav);
  set(gca,'ylim',[10.4 12.6],'xlim',[0 180],'ytick',10:12,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rates_cue1.eps']);
  
  % Plot DAN firing rate difference: RPE
  h = plot(dap-dav,'linewidth',0.2); set(h,'color',coldandiff);
  set(gca,'ylim',[-1.1 1.1],'xlim',[0 180],'ytick',-1:1,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rate_difference_cue1.eps']);
  close(gcf);
end;

if flag==3
  %%%%
  %%%% Fig. 2e & Supplementary Fig. 2
  %%%% VS model limited learning
  %%%% Dependence of learning on KC->DAN synaptic weights (gamma)
  %%%%
  
  %%% 5 values of gamma (w_k weights), each with 10 runs
  % Reward schedule
  rew = mb_reward_schedules(13,1,0,180,2);
  % Allocate memory
  m = zeros(180,10,5);
  go = zeros(180,10,5);
  nogo = zeros(180,10,5);
  dap = zeros(180,10,5);
  dav = zeros(180,10,5);
  wk = 0.8 + (0.1:0.1:0.5);
  % Generate data
  for j=1:10
    for k=1:5
      q = mb_vs(wk(k),j,45,13,2.5e-2,'rewardsign',1,'plasticity_rule',1,'choose1',true,'lambda',11.5);
      go(:,j,k) = q.go(:,1); nogo(:,j,k) = q.nogo(:,1);
      dap(:,j,k) = q.dap; dav(:,j,k) = q.dav;
      m(:,j,k) = q.go(:,1) - q.nogo(:,1);
    end;
  end;
  
  % Plot rewards and reward predictions
  myfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  for k=1:5
    h = plot(m(:,:,k),'linewidth',0.2); set(h,'color',colv1/2^((k-1)/3));
  end;
  set(gca,'ylim',[-2.2 2.2],'xlim',[0 180],'ytick',0-2:2:2,'xtick',[0 20:20:180],'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  myfig(2,7.5);
  for k=[1 3 5]
    h = plot(go(:,:,k),'linewidth',0.2); set(h,'color',colm1);
    hold on;
    h = plot(nogo(:,:,k),'linewidth',0.2); set(h,'color',colm2);
    set(gca,'ylim',[-0.1 2.6],'xlim',[0 180],'ytick',0:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_cue1_wk_' num2str(wk(k)) '.eps']);
  end;
  close(gcf);
  
  % Plot DAN firing rates
  myfig(2,7.5);
  for k=[1 3 5]
    h = plot(dap(:,:,k),'linewidth',0.2); set(h,'color',coldap);
    hold on;
    h = plot(dav(:,:,k),'linewidth',0.2); set(h,'color',coldav);
    if k==1
      set(gca,'ylim',[9 13],'xlim',[0 180],'ytick',9:2:13,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    elseif k==3
      set(gca,'ylim',[10 14],'xlim',[0 180],'ytick',10:2:14,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    elseif k==5
      set(gca,'ylim',[12 16],'xlim',[0 180],'ytick',12:2:16,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    end;
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rates_cue1_wk_' num2str(wk(k)) '.eps']);
  end;
  
  % Plot DAN firing rate difference: RPE
  for k=[1 3 5]
    plot(dap(:,:,k)-dav(:,:,k),'linewidth',0.2,'color',coldandiff);
    hold on;
    set(gca,'ylim',[-2.3 2.3],'xlim',[0 180],'ytick',-2:2:2,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rate_difference_cue1_wk_' num2str(wk(k)) '.eps']);
  end;
  close(gcf);
end;

if flag==4
  %%%%
  %%%% Supplementary Fig. 3
  %%%% VS model with inhibitory rewards (illustrated in Fig. 3a)
  %%%%
  
  %%% Create 10 runs
  % Reward schedule
  rew = mb_reward_schedules(13,1,0,180,2);
  % Allocate memory
  m = zeros(180,10);
  go = zeros(180,10);
  nogo = zeros(180,10);
  dap = zeros(180,10);
  dav = zeros(180,10);
  % Generate data
  for j=1:10
    q = mb_vs(1,j,45,13,2.5e-2,0,0,1,'choose1',true); % Setting variable decm = 10.5
    go(:,j) = q.go(:,1); nogo(:,j) = q.nogo(:,1);
    dap(:,j) = q.dap; dav(:,j) = q.dav;
    m(:,j) = q.go(:,1) - q.nogo(:,1);
  end;
  
  % Plot rewards and reward predictions
  myfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  h = plot(m,'linewidth',0.2); set(h,'color',colv1);
  set(gca,'ylim',[-2.2 2.2],'xlim',[0 180],'ytick',0-2:2:2,'xtick',[0 20:20:180],'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  myfig(2,7.5);
  h = plot(go,'linewidth',0.2); set(h,'color',colm1);
  hold on;
  h = plot(nogo,'linewidth',0.2); set(h,'color',colm2);
  set(gca,'ylim',[0 2],'xlim',[0 180],'ytick',0:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_cue1.eps']);
  close(gcf);
  
  % Plot DAN firing rates
  myfig(2,7.5);
  h = plot(dap,'linewidth',0.2); set(h,'color',coldap);
  hold on;
  h = plot(dav,'linewidth',0.2); set(h,'color',coldav);
  set(gca,'ylim',[8.6 11.3],'xlim',[0 180],'ytick',[9 11],'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rates_cue1.eps']);
  
  % Plot DAN firing rate difference: RPE
  h = plot(dap-dav,'linewidth',0.2); set(h,'color',coldandiff);
  set(gca,'ylim',[-1.2 1.2],'xlim',[0 180],'ytick',-1:1,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rate_difference_cue1.eps']);
  close(gcf);
end;

if flag==5
  %%%%
  %%%% Fig. 3d & Supplementary Fig. 4
  %%%% MV model, changing KC->DAN synaptic weights (gamma)
  %%%%
  
  %%% Create 10 runs
  % Reward schedule
  rew = mb_reward_schedules(13,1,0,180,2);
  % Allocate memory
  m = zeros(180,10,5);
  go = zeros(180,10,5);
  nogo = zeros(180,10,5);
  dap = zeros(180,10,5);
  dav = zeros(180,10,5);
  wk = 0:0.3:1.5;
  % Generate data
  for j=1:10
    for k=1:5
      q = mb_mv_a(wk(k),j,45,13,0.5*2.5e-2,'choose1',true,'no',2);
      go(:,j,k) = q.go(:,1); nogo(:,j,k) = q.nogo(:,1);
      dap(:,j,k) = q.dap; dav(:,j,k) = q.dav;
      m(:,j,k) = q.go(:,1) - q.nogo(:,1);
    end;
  end;
  
  % Plot rewards and reward predictions
  myfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  for k=1:5
    h = plot(m(:,:,k),'linewidth',0.2); set(h,'color',colv1/2^((k-1)/3));
  end;
  set(gca,'ylim',[-2.2 2.2],'xlim',[0 180],'ytick',0-2:2:2,'xtick',[0 20:20:180],'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  % Plot rewards and reward predictions inset
  myfig(1.75,1.75);
  for k=1:5
    h = plot(m(:,:,k),'linewidth',0.2); set(h,'color',colv1/2^((k-1)/3));
    hold on;
  end;
  set(gca,'ylim',[1.1 2.1],'xlim',[41 60],'ytick',[],'xtick',[],'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions_zoom.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  myfig(2,7.5);
  for k=[1 3 5]
    h = plot(go(:,:,k),'linewidth',0.2); set(h,'color',colm1);
    hold on;
    h = plot(nogo(:,:,k),'linewidth',0.2); set(h,'color',colm2);
    set(gca,'ylim',[-0.1 2.6],'xlim',[0 180],'ytick',0:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_cue1_wk_' num2str(wk(k)) '.eps']);
  end;
  close(gcf);
  
  % Plot DAN firing rates
  myfig(2,7.5);
  for k=[1 3 5]
    h = plot(dap(:,:,k),'linewidth',0.2); set(h,'color',coldap);
    hold on;
    h = plot(dav(:,:,k),'linewidth',0.2); set(h,'color',coldav);
    if k==1
      set(gca,'ylim',[0 1.5],'xlim',[0 180],'ytick',0:0.5:1.5,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    elseif k==3
      set(gca,'ylim',[4.5 7.5],'xlim',[0 180],'ytick',5:7,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    elseif k==5
      set(gca,'ylim',[10.5 13.5],'xlim',[0 180],'ytick',11:13,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    end;
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rates_cue1_wk_' num2str(wk(k)) '.eps']);
  end;
  close(gcf);
  
  % Plot DAN firing rate difference: RPE
  myfig(2,7.5);
  for k=[1 3 5]
    plot(dap(:,:,k)-dav(:,:,k),'linewidth',0.2,'color',coldandiff);
    hold on;
    set(gca,'ylim',[-2.6 2.6],'xlim',[0 180],'ytick',-2:2:2,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rate_difference_cue1_wk_' num2str(wk(k)) '.eps']);
  end;
  close(gcf);
end;

if flag==6
  %%%%
  %%%% Fig. 3d-g & Supplementary Fig. 4
  %%%% MV model
  %%%%
  
  % Create 10 runs
  % Reward schedule
  rew = mb_reward_schedules(13,1,0,180,2);
  % Allocate memory
  m = zeros(180,10);
  go = zeros(180,10);
  nogo = zeros(180,10);
  dap = zeros(180,10);
  dav = zeros(180,10);
  % Generate data
  for j=1:10
    q = mb_mv_a(1,j,45,13,0.5*2.5e-2,'choose1',true); % Setting variable decm = 10.5
    go(:,j) = q.go(:,1); nogo(:,j) = q.nogo(:,1);
    dap(:,j) = q.dap; dav(:,j) = q.dav;
    m(:,j) = q.go(:,1) - q.nogo(:,1);
  end;
  
  % Plot rewards and reward predictions
  myfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  h = plot(m,'linewidth',0.2); set(h,'color',colv1);
  set(gca,'ylim',[-2.4 2.4],'xlim',[0 180],'ytick',-2:2:2,'xtick',0:20:180,'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  qqmyfig(2,7.5);
  h = plot(go,'linewidth',0.2); set(h,'color',colm1);
  hold on;
  h = plot(nogo,'linewidth',0.2); set(h,'color',colm2);
  set(gca,'ylim',[0 2.2],'xlim',[0 180],'ytick',0:2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_cue1.eps']);
  close(gcf);
  
  % Plot DAN firing rates
  myfig(2,7.5);
  h = plot(dap,'linewidth',0.2); set(h,'color',coldap);
  hold on;
  h = plot(dav,'linewidth',0.2); set(h,'color',coldav);
  set(gca,'ylim',[8.6 11.3],'xlim',[0 180],'ytick',[9 11],'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rates_cue1.eps']);
  % Plot DAN firing rate difference: RPE
  h = plot(dap-dav,'linewidth',0.2); set(h,'color',coldandiff);
  set(gca,'ylim',[-2.6 2.6],'xlim',[0 180],'ytick',-2:2:2,'xtick',0:20:180,'xticklabel',{'0' [] [] '60' [] [] '120' [] [] '180'},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'DAN_firing_rate_difference_cue1.eps']);
  close(gcf);
end;

%%%
%%% Compare MV and VS-lambda (limited learning) models when
%%% each neuron class is blocked (Shibire) or activated (dTrpA1) during
%%% different stages of the simulation, for both appetitive (+ve reward)
%%% and aversive (-ve reward) training paradigms. 4 intervention protocols
%%% x 4 neuron types x 2 genetic interventions x 3 reward types in total:
%%% - Intervene during: CS+ only/CS+ & CS-/test only/throughout training & test
%%% - Target neuron: M+/M-/D+/D-
%%% - Genetic intervention: Shibire (block)/dTrpA1(activate)
%%% - Reward type: aversive/appetitive/none

if flag==8
  %%%%
  %%%% Generate associative conditioning experiment data for MV model
  %%%%
  
  % # RNG seeds
  nrep = 1000;
 
  % Setup parallel pool
   if oldflag    
     p = parcluster;
     NC = min(nrep,p.NumWorkers - 1);
     p = matlabpool('size');
     if p<1
       matlabpool(NC);
     else
       matlabpool close force local;
       matlabpool(NC);
     end;
  else
    p = gcp('nocreate');
    NC = min(nrep,7);
    if isempty(p)
      pp = parpool(NC);
    else
      if p.NumWorkers~=NC
        delete(p);
        pp = parpool(NC);
      end;
    end;
  end;
  
  % Init params
  nip = 4; % # intervention protocols (intervene during CS+ only/CS+ and CS-/testing only/CS+, CS- and testing)
  ntarg = 4; % # target neurons
  ngi = 2; % # genetic interventions (Shibire/TrpA)
  nrt = 3; % # reward types (+ve/-ve/none)
  nt = 30; % # trials
  nav = 50; % # runs over each block of averages
  nbl = nrep / nav; % # blocks of averages
  nit = floor(nrep/NC) + 1; % # iterations of spmd
  
  intrvn_id = [1 2 3 4];
  intrvn_amp = [0.1 5];
  rew_amp = [-1 1 0];
  
  % Allocate memory for interventions
  pre = zeros(nip,ntarg,ngi,nrt,nrep,nt,2); % Value predictions (m+ - m-)
  dec = zeros(nip,ntarg,ngi,nrt,nrep,nt); % Decisions
  % Allocate memory for controls
  cpre = zeros(nrt,nrep,nt,2); % Value predictions (m+ - m-)
  cdec = zeros(nrt,nrep,nt); % Decisions
  
  tic;
  % Run interventions
  for rt=1:nrt
    for ip=1:nip
      for targ=1:ntarg
        for gi=1:ngi
          for it=1:nit
            spmd(min(NC,nrep-(it-1)*NC))
              if gi==1
                intrvn_schedule = ones(nt,1);
              elseif gi==2
                intrvn_schedule = zeros(nt,1);
              end;
              % Generate intervention schedules for all trials
              if ip==1 % Intervention during CS+ training only
                intrvn_schedule(1:nt/3) = intrvn_amp(gi) * ones(nt/3,1);
              elseif ip==2 % Intervention during both CS+ training and CS-
                intrvn_schedule(1:2*nt/3) = intrvn_amp(gi) * ones(2*nt/3,1);
              elseif ip==3 % Intervention during test only
                intrvn_schedule((2*nt/3 + 1):nt) = intrvn_amp(gi) * ones(nt/3,1);
              elseif ip==4 % Intervention during both CS+, CS- and testing
                intrvn_schedule = intrvn_amp(gi) * ones(nt,1);
              end;
              % Generate reward schedule
              rew = mb_reward_schedules(12,(it-1)*NC+labindex,0.1,nt,2,rew_amp(rt));
              % Run model
              q = mb_mv_conditioning((it-1)*NC+labindex,rew,1.25e-2,'intervene_id',intrvn_id(targ),mod(gi,2),intrvn_schedule);
            end;
            for j=1:min(NC,nrep-(it-1)*NC)
              qq = q{j};
              pre(ip,targ,gi,rt,(it-1)*NC+j,:,1) = qq.map(:,1) - qq.mav(:,1);
              pre(ip,targ,gi,rt,(it-1)*NC+j,:,2) = qq.map(:,2) - qq.mav(:,2);
              dec(ip,targ,gi,rt,(it-1)*NC+j,:) = qq.decision;
            end;
            clear intrvn_schedule rew q;
          end;
          fprintf('rt: %d/%d,     ip: %d/%d,     targ: %d/%d,     gi: %d/%d\n',rt,nrt,ip,nip,targ,ntarg,gi,ngi);
        end;
      end;
    end;
  end;
  toc
  
  % Run controls
  tic;
  for rt=1:2 % Don't run controls for experiments with no reinforcement
    for it=1:nit
      spmd(min(NC,nrep-(it-1)*NC))
        % Generate reward schedule
        rew = mb_reward_schedules(12,(it-1)*NC+labindex,0.1,nt,2,rew_amp(rt));
        % Run model
        q = mb_mv_conditioning((it-1)*NC+labindex,rew,0.5*2.5e-2);
      end;
      for j=1:min(NC,nrep-(it-1)*NC)
        qq = q{j};
        cpre(rt,(it-1)*NC+j,:,1) = qq.map(:,1) - qq.mav(:,1);
        cpre(rt,(it-1)*NC+j,:,2) = qq.map(:,2) - qq.mav(:,2);
        cdec(rt,(it-1)*NC+j,:) = qq.decision;
      end;
      clear rew q;
    end;
  end;
  toc
  save([fpath 'associative_conditioning_data_for_MV_model.mat'],'pre','dec','cpre','cdec');
end;

if flag==9
  %%%%
  %%%% Generate associative conditioning experiment data for VS-lambda model
  %%%%
  
  % # RNG seeds
  nrep = 1000;
  
   % Setup parallel pool
   if oldflag    
     p = parcluster;
     NC = min(nrep,p.NumWorkers - 1);
     p = matlabpool('size');
     if p<1
       matlabpool(NC);
     else
       matlabpool close force local;
       matlabpool(NC);
     end;
  else
    p = gcp('nocreate');
    NC = min(nrep,7);
    if isempty(p)
      pp = parpool(NC);
    else
      if p.NumWorkers~=NC
        delete(p);
        pp = parpool(NC);
      end;
    end;
  end;

  % Init params
  nip = 4; % # intervention protocols (intervene during training/testing/both)
  ntarg = 4; % # target neurons
  ngi = 2; % # genetic interventions (Shibire/TrpA)
  nrt = 3; % # reward types (+ve/-ve/none)
  nt = 30; % # trials
  nav = 50; % # runs over each block of averages
  nbl = nrep / nav; % # blocks of averages
  nit = floor(nrep/NC) + 1; % # iterations of spmd
  
  intrvn_id = [1 2 3 4];
  intrvn_amp = [0.1 5];
  rew_amp = [-1 1 0];
  
  % Allocate memory for interventions
  pre = zeros(nip,ntarg,ngi,nrt,nrep,nt,2); % Value predictions (m+ - m-)
  dec = zeros(nip,ntarg,ngi,nrt,nrep,nt); % Decisions
  % Allocate memory for controls
  cpre = zeros(nrt,nrep,nt,2); % Value predictions (m+ - m-)
  cdec = zeros(nrt,nrep,nt); % Decisions
  
  tic;
  % Run interventions
  for rt=1:nrt
    for ip=1:nip
      for targ=1:ntarg
        for gi=1:ngi
          for it=1:nit
            spmd(min(NC,nrep-(it-1)*NC))
              if gi==1
                intrvn_schedule = ones(nt,1);
              elseif gi==2
                intrvn_schedule = zeros(nt,1);
              end;
              if ip==1 % Intervention during CS+ training only
                intrvn_schedule(1:nt/3) = intrvn_amp(gi) * ones(nt/3,1);
              elseif ip==2 % Intervention during both CS+ training and CS-
                intrvn_schedule(1:2*nt/3) = intrvn_amp(gi) * ones(2*nt/3,1);
              elseif ip==3 % Intervention during test only
                intrvn_schedule((2*nt/3 + 1):nt) = intrvn_amp(gi) * ones(nt/3,1);
              elseif ip==4 % Intervention during both CS+, CS- and testing
                intrvn_schedule = intrvn_amp(gi) * ones(nt,1);
              end;
              rew = mb_reward_schedules(12,(it-1)*NC+labindex,0.1,nt,2,rew_amp(rt));
              q = mb_vs_conditioning((it-1)*NC+labindex,rew,5e-2,'intervene_id',intrvn_id(targ),mod(gi,2),intrvn_schedule);
            end;
            for j=1:min(NC,nrep-(it-1)*NC)
              qq = q{j};
              pre(ip,targ,gi,rt,(it-1)*NC+j,:,1) = qq.map(:,1) - qq.mav(:,1);
              pre(ip,targ,gi,rt,(it-1)*NC+j,:,2) = qq.map(:,2) - qq.mav(:,2);
              dec(ip,targ,gi,rt,(it-1)*NC+j,:) = qq.decision;
            end;
            clear intrvn_schedule rew q;
          end;
          fprintf('rt: %d/%d,     ip: %d/%d,     targ: %d/%d,     gi: %d/%d\n',rt,nrt,ip,nip,targ,ntarg,gi,ngi);
        end;
      end;
    end;
  end;
  toc
  
  % Run controls
  tic;
  for rt=1:2 % Don't run controls for experiments with no reinforcement
    for it=1:nit
      spmd(min(NC,nrep-(it-1)*NC))
        % >>> reward schedule needs rng seed, ID should =12
        rew = mb_reward_schedules(12,(it-1)*NC+labindex,0.1,nt,2,rew_amp(rt));
        q = mb_vs_conditioning((it-1)*NC+labindex,rew,5e-2);
      end;
      for j=1:min(NC,nrep-(it-1)*NC)
        qq = q{j};
        cpre(rt,(it-1)*NC+j,:,1) = qq.map(:,1) - qq.mav(:,1);
        cpre(rt,(it-1)*NC+j,:,2) = qq.map(:,2) - qq.mav(:,2);
        cdec(rt,(it-1)*NC+j,:) = qq.decision;
      end;
      clear rew q;
    end;
  end;
  toc
  
  save([fpath 'associative_conditioning_data_for_VS_model.mat'],'pre','dec','cpre','cdec');
end;

if flag==10
  %%%%
  %%%% Fig. 5d-e, analyse VS_lambda and MV model learned behaviours in the
  %%%% 96 possible interventions plus controls in an associative
  %%%% conditioning paradigm. Load experimental data for comparison.
  %%%%
  
  m = cell(2,1);
  m{1} = load([fpath 'associative_conditioning_data_for_MV_model.mat']);
  m{2} = load([fpath 'associative_conditioning_data_for_VS_model.mat']);
  nrep = 1000; % # repeat runs of the simulation per condition
  nip = 4; % # intervention protocols (intervene during training/testing/both)
  ntarg = 4; % # target neurons
  ngi = 2; % # genetic interventions (Shibire/TrpA)
  nrt = 3; % # reward types (+ve/-ve)
  nt = 30; % # trials
  nav = 50; % # runs over each block of averages
  nbl = nrep / nav; % # blocks of averages
  nmod = 3; % # models
  cc = cell(2,1); % Cell array containg performance index (PI) data from simulated genetic intervention conditions
	ct = cell(2,1); % Cell array containg performance index (PI) data from simulated control conditions
  plmod = cell(2,1); % Cell array for plotting data from the model
  plexp = cell(2,1); % Cell array for plotting data from experiments
  
  %%% Orgainse model PI data (experimental conditions) for further analysis
  for model=1:2
    condpi = zeros(nip,ntarg,ngi,nrt,nbl,nt/2);
    for ip=1:nip
      for targ=1:ntarg
        for gi=1:ngi
          for bl=1:nbl
            for rt=1:nrt
              for t=1:nt/2
                temp = squeeze(m{model}.dec(ip,targ,gi,rt,(bl-1)*nav+1:bl*nav,(t-1)*2+1:t*2));
                n1 = sum(temp(:)==1); n2 = sum(temp(:)==2);
                condpi(ip,targ,gi,rt,bl,t) = (n1-n2) / (n1+n2);
              end;end;end;end;end;end;
    ccondpi = reshape(condpi,nip*ntarg*ngi*nbl*nrt,nt/2);
    ccondpi(:,1:10) = [];
    cc{model} = ccondpi;
  end;
  
  % Organise model PI data (controls) for further analysis
	for model=1:2
		contpi = zeros(2,nbl,2,nt/2);
		for rt=1:2
			for bl=1:nbl
				for j=1:nt/2
					temp = squeeze(m{model}.cdec(rt,(bl-1)*nav+1:bl*nav,(j-1)*2+1:j*2));
					n1 = sum(temp(:)==1); n2 = sum(temp(:)==2);
					contpi(rt,bl,j) = (n1-n2) / (n1+n2);
				end;end;end;
		ccontpi = reshape(contpi,rt*nbl,nt);
		ccontpi(:,1:10) = [];
		ct{model} = ccontpi;
	end;
  
  % Generate indeces into conditions in cell array 'cc'
  a = zeros(nip,ntarg,ngi,nrt); j = 1;
  for ip=1:nip
    for targ=1:ntarg
      for gi=1:ngi
        for rt=1:nrt
          a(ip,targ,gi,rt) = j;
          j = j + 1;
        end;end;end;end;
  a = repmat(a,[1,1,1,1,nbl,1]);
  indcond = reshape(a,nip*ntarg*ngi*nrt*nbl,1);
  
  % Generate a list of all possible conditions, where each entry in the
  % list is a 1x4 vector (hereafter called a condition key). In each
  % vector, the 1st element denotes the
  % intervention protocol, the 2nd element denotes the target neuron
  % type, 3rd element denotes the genetic intervention (shi/dTrpA1), the
  % 4th element denotes the reward type.
  condslist = zeros(nip*ntarg*ngi*nrt,4);
  j=1;
  for ip=1:nip
    for targ=1:ntarg
      for gi=1:ngi
        for rt=1:nrt
          condslist(j,:) = [ip targ gi rt];
          j = j + 1;
        end;end;end;end;
  
  % Read in experimental PI measures from Excel sheet
	z = xlsread([fpath 'supptable2.xlsx'],'A3:AG167');

  testtimes = z(:,end);
  z = z(:,[1 30]);
  % Compute indices into the rows of the experimental data array where
  % there is data
  ind0 = find(~isnan(z(:,2)));
  %     % Compute indices into the rows of the experimental data array where
  %     % there is data and where the test stage was performed within a desired
  %     % time limit.
  %     ind0 = find(~isnan(z(:,2))&(testtimes<31));
  
  % Generate a list of condition keys for each row of the experimental
  % data array (this may contain repeats of each key listed in
  % 'condslist').
  ex = zeros(length(ind0),4);
  for j=1:length(ind0)
    scond = num2str(z(ind0(j),1));
    for k=1:4
      ex(j,k) = str2num(scond(k));
    end;
  end;
  
  % Generate different indices into rows of experimental data array,
  % according to specific cases of each of the 4 experimental parameters.
  ind0ip = false(length(ind0),4); % Indices for each of the 4 intervention protocols
  ind0targ = false(length(ind0),4); % Indices for each of the 4 target neurons
  ind0gi = false(length(ind0),2); % Indices for each of the 2 genetic interventions
  ind0rt = false(length(ind0),3); % Indices for each of the 3 reward types
  for j=1:4
    ind0ip(:,j) = ex(:,1)==j;
    ind0targ(:,j) = ex(:,2)==j;
  end;
  for j=1:2
    ind0gi(:,j) = ex(:,3)==j;
  end;
  for j=1:3
    ind0rt(:,j) = ex(:,4)==j;
  end;
  
  % Populate the cell arrays for plotting
  for model=1:2
    plmod{model} = zeros(nbl,length(ind0));
    plexp{model} = zeros(nbl,length(ind0));
    for j=1:length(ind0)
      if ex(j,4)==1
        temp = ct{model}(1:2:end,1);
        cont = mean(temp(:))*ones(nbl,1);
      elseif ex(j,4)==2
        temp = ct{model}(2:2:end,1);
        cont = mean(temp(:))*ones(nbl,1);
      elseif ex(j,4)==3
        cont = zeros(nbl,1);ones(nbl,1);
      end;
      
      temp = (cc{model}(indcond==find(ismember(condslist,ex(j,:),'rows')),1) + 1) / 2;
      cont = (cont + 1) / 2;
      plmod{model}(:,j) = (temp - cont) ./ sqrt((temp + cont)/2.*(1 - (temp+cont)/2)*2/nav);
      indnan = isnan(plmod{model}(:,j));
      plmod{model}(indnan,j) = 0;
      % Make repeat copies of each experimetnal data point, to pair up
      % with the 20 model data points for each condition
      plexp{model}(:,j) = z(ind0(j),2) * ones(nbl,1);
    end;
  end;
  
  % Generate cell array of colours for plotting
  colneurons = cell(4,1);
  for j=1:4
    switch true
      case j==1
        colneurons{j} = colm1;
      case j==2
        colneurons{j} = colm2;
      case j==3
        colneurons{j} = coldap;
      case j==4
        colneurons{j} = coldav;
    end;
  end;
  
  % Plot the experimental versus model test statistics for differences
  % between control PIs and condition PIs.  
  for j=1:2
    myfig(8,8);    
    % Compute the weighted linear fit
    [b s] = robustfit(plmod{j}(:),plexp{j}(:));
    % Weights from robustfit
    weights = reshape(s.w,nbl,length(ind0));
    % Plot, coloured by neuron type    
    for m=[1 3 4 2]
      for n=1:2
        ind = find(ind0targ(:,m)&ind0gi(:,n));
        for k=1:length(ind)
          for q=1:nbl
            plot(plmod{j}(q,ind(k)),plexp{j}(q,ind(k))+0.5*(rand-.5),'o','color',colneurons{m},'markersize',1+2*weights(q,ind(k)));
            hold on;
          end;
        end;
      end;
    end;    
    % Plot weighted linear fit
    x = [-6.6 6.6]; plot(x,b'*[1,1;-6.6,6.6],'linewidth',1,'color',[0.5 0.5 0.5]);
    set(gca,'ylim',[-14 10],'xlim',[-6.6 6.6],'ytick',-10:5:10,'xtick',-6:3:6,'box','off');
    hold off;
    % Print out Pearson correlation coefficient
    fprintf(['Model' num2str(j) ':Corr coef Least Squares = %f\n'],corr(plmod{j}(:),plexp{j}(:)));
    fprintf(['Model' num2str(j) ':Corr coef Robust Fit = %f\n'],corr(plmod{j}(:).*s.w,plexp{j}(:).*s.w));
    myticks;
    if j==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rho_model_dap_dav_vs_exp_pi_diff.eps']);
    elseif j==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rho_model_decm_da_vs_exp_pi_diff.eps']);
    end;
    delete(gcf);
  end;
  
  % Plot examples from above that will be highlighted in the text
  exampleind = [33 53 85 89, 1 8 22]; 
  for j=1:2
    myfig(8,8);    
    % Compute the weighted linear fit
    [b s] = robustfit(plmod{j}(:),plexp{j}(:));
    % Weights from robustfit
    weights = reshape(s.w,nbl,length(ind0));
    % Plot first set of examples, coloured by neuron type.
    for k=1:4
      ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
      ind = find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi));
      for m=1:length(ind)
        for q=1:nbl
          plot(plmod{j}(q,ind(m)),plexp{j}(q,ind(m))+0.2*(rand-.5),'o','color',colneurons{k},'markersize',1+2*weights(q,ind(m)));
          hold on;
        end;
      end;
    end;hold off;    
    set(gca,'ylim',[-14 10],'xlim',[-6.6 6.6],'ytick',-10:5:10,'xtick',-6:3:6,'box','off');
    myticks;
    if j==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rho_model_dap_dav_vs_exp_examples.eps']);
    elseif j==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rho_model_decm_da_vs_exp_examples.eps']);
    end;
    % Plot second set of examples, coloured by neuron type.
    for k=5:7
      ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
      ind = find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi));
      for q=1:nbl
        for m=1:length(ind)
          plot(plmod{j}(q,ind(m)),plexp{j}(q,ind(m))+0.2*(rand-.5),'o','color',(colneurons{targ,1}-0.05).^2,'markersize',1+2*weights(q,ind(m)));
          hold on;
        end;
      end;
    end; hold off;
    set(gca,'ylim',[-14 10],'xlim',[-6.6 6.6],'ytick',-10:5:10,'xtick',-6:3:6,'box','off');
    myticks;
    if j==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rho_model_dap_dav_vs_exp_examples2.eps']);
    elseif j==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rho_model_decm_da_vs_exp_examples2.eps']);
    end;
    delete(gcf);
  end;
end;

if flag==11
  %%%%
  %%%% Fig. 4a, example of the MV model's performance for 10 runs of a
  %%%% 2-alternative forced choice (2-armed bandit) task
  %%%%
  % Reward schedule
  rew = mb_reward_schedules(3,20,0,250,2,2,10); rew(1:50,:) = [];
  % Generat4e data for 10 runs, using same reward schedule
  pred = zeros(200,2,10);
  dec = zeros(200,10);
  for j=1:10
    q = mb_mv_a(1,j,200,rew,2.5e-2,'no',2);
    pred(:,:,j) = q.go - q.nogo;
    dec(:,j) = q.decision;
  end;
  
  % Plot rewards
  myfig(4,7.5);
  plot(rew(:,1),'linewidth',2,'color',colr1);
  hold on;
  plot(rew(:,2),'linewidth',2,'color',colr2);  
  set(gca,'ylim',[-3 2.1],'xlim',[0 200],'ytick',-2:2:2,'xtick',0:100:200,'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards.eps']);
  % Plot reward predictions
  for j=1:10
    plot(find(dec(:,j)==1),pred(dec(:,j)==1,1,j),'o','color',colv1,'markersize',2,'linewidth',0.5);
    hold on;
    plot(find(dec(:,j)==2),pred(dec(:,j)==2,2,j),'o','color',colv2,'markersize',2,'linewidth',0.5);
  end; hold off;
  set(gca,'ylim',[-3 2.1],'xlim',[0 200],'ytick',-2:2:2,'xtick',0:100:200,'box','off');
  axis off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'reward_predictions.eps']);
  close(gcf);  
end;

if flag==12
  %%%%
  %%%% Fig. 4b, example of the MV model's performance for a
  %%%% 5-alternative forced choice (5-armed bandit) task
  %%%%
  
  %%% Generate data for a single run
  % Reward schedule  
  rew = mb_reward_schedules(3,20,0,250,5,2,10); rew(1:50,:) = [];
  % Run model
  q = mb_mv_a(1,20,200,rew,2.5e-2,'no',5);
  
  % Plot rewards
  myfig(4,7.5);
  for j=1:5
    col = ([((j-1)/4) abs(j-3)/3 ((5-j)/4)]);
    plot(q.r(:,j),'linewidth',2,'color',col);
    hold on;
  end; hold off;
  set(gca,'ylim',[-3 2.1],'xlim',[0 200],'ytick',-2:2:2,'xtick',0:100:200,'box','off');
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions_1.eps']);  
  % Plot reward predictions (chosen cues)
  for j=1:5
    col = ([((j-1)/4) abs(j-3)/3 ((5-j)/4)]);
    plot(find(q.decision==j),q.go(q.decision==j,j) - q.nogo(q.decision==j,j),'o','color',col,'markersize',2,'linewidth',0.5);
    hold on;
  end;hold off;
  set(gca,'ylim',[-3 2.1],'xlim',[0 200],'ytick',-2:2:2,'xtick',0:100:200,'box','off');
  axis off     
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions_2.eps']);
  close(gcf);
end;

if flag==13
  %%%%
  %%%% Supplementary Fig. 5a, testing VS_lambda and MV model performance 
  %%%% as a function of the number of cues (from 2-200).
  %%%%
  
  % Setup parallel pool
  nrep = 100;
   if oldflag    
     p = parcluster;
     NC = min(nrep,p.NumWorkers - 1);
     p = matlabpool('size');
     if p<1
       matlabpool(NC);
     else
       matlabpool close force local;
       matlabpool(NC);
     end;
  else
    p = gcp('nocreate');
    NC = min(nrep,7);
    if isempty(p)
      pp = parpool(NC);
    else
      if p.NumWorkers~=NC
        delete(p);
        pp = parpool(NC);
      end;
    end;
  end;
  
  % Initialise parameters
  nstim = [2 3 4 6 8 12 16 24 32 48 64 96 128 200]; % # stimuli to test
  % Allocate memory
  rewpercue1 = zeros(length(nstim),nrep); % For 100% accurate agent, beta = inf (black solid line)
  rewpercue2 = zeros(length(nstim),nrep); % For 100% accurate agent, beta = 5 (black dotted line)
  rewpercue3 = zeros(length(nstim),nrep); % For random agent (black dashed line)
  rewpercue4 = cell(4,1); % Perfect plasticity agent and MV model
  for j=1:4
    rewpercue4{j} = zeros(length(nstim),nrep);
  end;
  nit = floor(nrep/NC) + 1; % # iterations of spmd
  tic;
  for model=1:4
    for no=1:length(nstim)
      for it=1:nit
        spmd(min(NC,nrep-(it-1)*NC))
          % Reward schedule
          rew = mb_reward_schedules(3,(it-1)*NC+labindex,0,250,nstim(no),2,10);
          rew(1:50,:) = [];
          
          mrew = sum(exp(rew/0.2) ./ repmat(sum(exp(rew / 0.2),2),1,nstim(no)) .* rew,2); % Mean reward (reward x probabiity of choosing rewarding cue (from softmax))
          if model==1
            % Overlapping KC representations, MV model
            q = mb_mv_a(1,(it-1)*NC+labindex,200,rew,2.5e-2,'no',nstim(no));
          elseif model==2
            % Non-overlapping KC representations, MV model
            q = mb_mv_d(1,(it-1)*NC+labindex,200,rew,2.5e-2,'no',nstim(no));
          elseif model==3
            % Non-overlapping KC representations, perfect synaptic weight updates
            q = mb_mv_c(1,(it-1)*NC+labindex,200,rew,'no',nstim(no));
          elseif model==4
            % Non-overlapping KC representations, VS model
%             q = mb_vs(1,(it-1)*NC+labindex,200,rew,2.5e-1,'no',nstim(no),'lambda',12,'nk',no*100);
            q = mb_vs(1,(it-1)*NC+labindex,200,rew,1e-1,'no',nstim(no),'lambda',12,'nk',nstim(no)*100);
          end;
        end;
        for j=1:min(NC,nrep-(it-1)*NC)
          qq = q{j};
          if model==1
            rewpercue1(no,(it-1)*NC+j) = sum(max(qq.r,[],2)) / 200; % Max reward per trial
            rewpercue2(no,(it-1)*NC+j) = sum(mrew{j}) / 200; % Weighted mean reward per trial
            rewpercue3(no,(it-1)*NC+j) = sum(mean(qq.r,2)) / 200; % Mean reward per trial
          end;
          obr = qq.sr; % Obtained reward per trial          
          rewpercue4{model}(no,(it-1)*NC+j) = obr / 200;
        end;
        clear rew q mrew;
        fprintf('model = %d/4,    no = %d/%d,    it = %d/%d\n',model,no,length(nstim),it,nit);
      end;
    end;    
  end;
  toc
  % Plot SDs and means for the different agents and MV models
  myfig(4,7.5);  
  col = [191, 84, 64; 26, 128, 204; 105, 0, 190; 128, 204, 26] / 255;
  for model=1:4
    if model==1
      h = fill([nstim'; fliplr(nstim)'],[mean(rewpercue1,2)+std(rewpercue1,1,2); flipud(mean(rewpercue1,2)-std(rewpercue1,1,2))],'k');
      set(h,'edgecolor','none','facealpha',0.35);
      hold on;
      plot(nstim',mean(rewpercue1,2),'k','linewidth',1);
      h = fill([nstim'; fliplr(nstim)'],[mean(rewpercue2,2)+std(rewpercue2,1,2); flipud(mean(rewpercue2,2)-std(rewpercue2,1,2))],'k');
      set(h,'edgecolor','none','facealpha',0.35);
      plot(nstim',mean(rewpercue2,2),'linewidth',1,'color','k');
      h = fill([nstim'; fliplr(nstim)'],[mean(rewpercue3,2)+std(rewpercue3,1,2); flipud(mean(rewpercue3,2)-std(rewpercue3,1,2))],'k');
      set(h,'edgecolor','none','facealpha',0.35);
      plot(nstim',mean(rewpercue3,2),'linewidth',1,'color','k');
      h1 = fill([nstim'; fliplr(nstim)'],[mean(rewpercue4{model},2)+std(rewpercue4{model},1,2); flipud(mean(rewpercue4{model},2)-std(rewpercue4{model},1,2))],col(model,:));
      set(h1,'edgecolor','none','facealpha',0.35); hold on;
      h2 = plot(nstim',mean(rewpercue4{model},2),'linewidth',1,'color',col(model,:));
      set(gca,'xscale','log','ylim',[-0.35 1.7],'ytick',0:0.5:1.5,'xlim',[2 200],'xtick',[2:2:10 20:20:100 200],'xticklabel',{2 [] [] [] 10 20 [] [] [] 100 200},'box','off');
      myticks; hold off;
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'reward_per_trial_overlappingcues_MV_model.eps']);      
    else
      h1 = fill([nstim'; fliplr(nstim)'],[mean(rewpercue4{model},2)+std(rewpercue4{model},1,2); flipud(mean(rewpercue4{model},2)-std(rewpercue4{model},1,2))],col(model,:));
      set(h1,'edgecolor','none','facealpha',0.35); hold on;
      h2 = plot(nstim',mean(rewpercue4{model},2),'linewidth',1,'color',col(model,:));
      set(gca,'xscale','log','ylim',[-0.35 1.7],'ytick',0:0.5:1.5,'xlim',[2 200],'xtick',[2:2:10 20:20:100 200],'xticklabel',{2 [] [] [] 10 20 [] [] [] 100 200},'box','off');
      myticks; hold off;
      if model==2
        myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'reward_per_trial_nonoverlappingcues_MV_model.eps']);
      elseif model==3
        myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'reward_per_trial_perfect_pred_nonoverlappingcues.eps']);
      elseif model==4
        myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'reward_per_trial_nonoverlappingcues_VS_model.eps']);
      end;
    end;
    delete(h1); delete(h2);
  end;
  delete(gcf);
end;

if flag==14
  %%%%
  %%%% Supplementary Fig. 5b, testing differences in performance between
  %%%% the perfect plasticity agent and the MV model as a function of
  %%%% noise in the rewards and the # of cues.
  %%%%
  
  nrep = 100;
  % Setup parallel pool
   if oldflag    
     p = parcluster;
     NC = min(nrep,p.NumWorkers - 1);
     p = matlabpool('size');
     if p<1
       matlabpool(NC);
     else
       matlabpool close force local;
       matlabpool(NC);
     end;
  else
    p = gcp('nocreate');
    NC = min(nrep,7);
    if isempty(p)
      pp = parpool(NC);
    else
      if p.NumWorkers~=NC
        delete(p);
        pp = parpool(NC);
      end;
    end;
  end;
  
  nit = floor(nrep/NC) + 1; % # iterations of spmd
  
  nstim = [4 20 40]; % # cues
  obr = zeros(nrep,11,2,length(nstim));
  for no=1:3
    for noise=1:11
      for it=1:nit
        spmd(min(NC,nrep-(it-1)*NC))
          rew = mb_reward_schedules(3,(it-1)*NC+labindex,(noise-1)*0.3,250,nstim(no),2,10);
          rew(1:50,:) = [];
          q1 = mb_mv_d(1,(it-1)*NC+labindex,200,rew,5e-2,'no',nstim(no));
          q2 = mb_mv_c(1,(it-1)*NC+labindex,200,rew,'no',nstim(no));
        end;
        for j=1:min(NC,nrep-(it-1)*NC)
          qq1 = q1{j};
          qq2 = q2{j};
          obr((it-1)*NC+j,noise,1,no) = sum(sum(qq1.r(repmat(qq1.decision,1,nstim(no)) == ones(length(qq1.decision),1)*(1:nstim(no))))) / 200;
          obr((it-1)*NC+j,noise,2,no) = sum(sum(qq2.r(repmat(qq2.decision,1,nstim(no)) == ones(length(qq2.decision),1)*(1:nstim(no))))) / 200;
        end;
        clear rew q1 q2;
      end;
    end;
  end;
  
  myfig(4,7.5);
  col = [26, 128, 204; 105, 0, 190] / 255; % Coloursfor the two curves
  x = 0:0.3:3; % Noise amplitudes for plotting
  for no=1:3
    h = fill([x'; fliplr(x)'],[mean(obr(:,:,1,no),1)'+std(obr(:,:,1,no),1,1)'; flipud(mean(obr(:,:,1,no),1)'-std(obr(:,:,1,no),1,1)')],col(1,:));
    set(h,'edgecolor','none','facealpha',0.35);
    hold on;
    plot(x,mean(obr(:,:,1,no),1),'color',col(1,:),'linewidth',1);
    h = fill([x'; fliplr(x)'],[mean(obr(:,:,2,no),1)'+std(obr(:,:,2,no),1,1)'; flipud(mean(obr(:,:,2,no),1)'-std(obr(:,:,2,no),1,1)')],col(2,:));
    set(h,'edgecolor','none','facealpha',0.35);
    plot(x,mean(obr(:,:,2,no),1),'color',col(2,:),'linewidth',1);
    set(gca,'xlim',[0 3],'xtick',0:3,'ylim',[-0.2 0.9],'ytick',0:0.3:0.9,'box','off')
    myticks; hold off;
    myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'reward_per_trial_versus_reward_noise_' num2str(nstim(no)) 'cues.eps']);
  end;
  save([fpath 'flag' num2str(flag) 'reward_per_trial_versus_reward_noise.mat'],'obr','nstim');
end;

if flag==15
  %%%%
  %%%% Supplementary Fig. 5d-f, Comparing reward predictions for single
  %%%% KCs and KC populations with actual rewards for mb_mv_b, where
  %%%% predictions for all cues are updated.
  %%%%
  myfig(7,7);
  
  nstim = 200; % # stimuli
  seed = 1;
  % Reward schedule
  rew = mb_reward_schedules(3,seed,0,250,nstim,2,10);
  rew(1:50,:) = [];
  % Run model
  q = mb_mv_b(1,seed,50,rew,1e-3,nstim);
  
  %%% Plot panel d
  %%% Each data point: reward prediction contribution from a single KC
  %%% for a single cue, versus the actual reward for that cue
  colcl = zeros(5,3);
  colcl(1,:) = [130 0 0] / 255;
  colcl(2,:) = [255 0 195] / 255;
  colcl(3,:) = [45 0 150] / 255;
  colcl(4,:) = [0 145 255] / 255;
  colcl(5,:) = [0 150 80] / 255;
  for kc=1:2000
    plot(q.r(end,q.s(kc,:)>0),q.s(kc,q.s(kc,:)>0)'*squeeze(q.wkmap(1,kc,end)-q.wkmav(1,kc,end))','ok','markersize',1.5,'color','k');
    hold on;
  end;
  for j=1:5
    kc = j+90;
    plot(q.r(end,q.s(kc,:)>0),q.s(kc,q.s(kc,:)>0)'*squeeze(q.wkmap(1,kc,end)-q.wkmav(1,kc,end))','ok','markersize',1.5,'markerfacecolor',colcl(j,:),'color',colcl(j,:));
    hold on;
  end;
  plot([-2 2],[-2 2]/100,'r','linewidth',1);
  set(gca,'xlim',[-1.7 1.7],'xtick',-1:1,'ylim',[-7.5 7.5]/100,'ytick',(-6:3:6)/100,'box','off');
  hold off;myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'prediction_vs_reward_per_kc_per_cue.eps']);
  
  %%% Plot panel e
  %%% Each data point: reward prediction contribution from a single KC
  %%% averaged over all cues, versus the actual mean reward for all cues
  %%% that elicit a response in that KC
  for kc=1:2000
    plot(mean(q.r(end,q.s(kc,:)>0)),mean(q.s(kc,q.s(kc,:)>0)'*squeeze(q.wkmap(1,kc,end)-q.wkmav(1,kc,end))'),'ok','markersize',2);
    hold on;
  end;
  plot([-2 2],[-2 2]/100*(0.05*nstim),'r','linewidth',1);
  set(gca,'xlim',[-1.7 1.7],'xtick',-1:1,'ylim',[-7.5 7.5]/100,'ytick',(-6:3:6)/100,'box','off');
  hold off;myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'prediction_vs_reward_per_kc_averaged_over_cues.eps']);
  
  %%% Plot panel f
  %%% Each data point: actual reward predictions versus acutal rewards
  %%% for each cue.
  rp = zeros(2000,nstim);
  for kc=1:2000
    rp(kc,:) = q.s(kc,:)*(q.wkmap(1,kc,end)-q.wkmav(1,kc,end));
  end;
  rp(rp==0) = nan;
  plot(q.r(end,:),nansum(rp,1),'ok','markersize',2);
  set(gca,'xlim',[-1.5 1.5],'xtick',-1:1,'ylim',[-1.5 1.5],'ytick',-1:1,'box','off');
  hold on; plot([-2 2],[-2 2],'r','linewidth',1); hold off; myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'prediction_vs_reward_per_cue.eps']);
  
  delete(gcf);
end;

