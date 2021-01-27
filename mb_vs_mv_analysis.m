function mb_vs_mv_analysis(flag)
%
% MB_VS_MV_ANALYSIS(FLAG)
%
% Analyse data for Bennett et al. 2019.
% Valence specific (VS) and mixed valence (MV) MB models are analysed.

%%%
%%% Output directory, where plots will be saved
%%%
fpath = '<directory name>';
% fpath = '/Users/jamesbennett/Documents/DPhil/doc/applications/sussex/manuscripts/2018bennett.nowotny_mb_resc_wag/submission/';
fpath = '~/sussex/manuscripts/2018bennett.nowotny_mb_resc_wag/revision3/';
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
    q = mb_vs(1,j,180,13,2.5e-2,'rewardsign',1,'plasticity_rule',0,'choose1',true);
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
  wap = zeros(180,10,10);
  wav = zeros(180,10,10);
  % Generate data
  for j=1:10
    q = mb_vs(1,j,180,13,2.5e-2,'rewardsign',1,'plasticity_rule',1,'choose1',true,'lambda',11.5);
    go(:,j) = q.go(:,1); nogo(:,j) = q.nogo(:,1);
    dap(:,j) = q.dap; dav(:,j) = q.dav;
    m(:,j) = q.go(:,1) - q.nogo(:,1);
    wap(:,:,j) = squeeze(q.wkmap(1,1:10,:))';
    wav(:,:,j) = squeeze(q.wkmav(1,1:10,:))';
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
  
  % Plot synaptic weights
  myfig(4,7.5);
  for j=1:10
    plot(wap(:,:,j),'linewidth',0.2,'color',colm1);
    hold on;
    plot(wav(:,:,j),'linewidth',0.2,'color',colm2);
  end;
  set(gca,'ylim',[0 0.25],'xlim',[0 180],'ytick',0:0.1:0.2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'weights_ap_and_av.eps']);
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
      q = mb_vs(wk(k),j,180,13,2.5e-2,'rewardsign',1,'plasticity_rule',1,'choose1',true,'lambda',11.5);
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
    q = mb_vs(1,j,180,13,2.5e-2,0,0,1,'choose1',true); % Setting variable decm = 10.5
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
  %%%% Fig. 3d & Supplementary Fig. 5
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
      q = mb_mv_a(wk(k),j,180,13,0.25*2.5e-2,'choose1',true,'no',2);
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
  %%%% Fig. 3d-g & Supplementary Fig. 5
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
  wap = zeros(180,10,10);
  wav = zeros(180,10,10);
  
  % Generate data
  for j=1:10
    q = mb_mv_a(1,j,180,13,0.5*2.5e-2,'choose1',true);
    go(:,j) = q.go(:,1); nogo(:,j) = q.nogo(:,1);
    dap(:,j) = q.dap; dav(:,j) = q.dav;
    m(:,j) = q.go(:,1) - q.nogo(:,1);
    wap(:,:,j) = squeeze(q.wkmap(1,1:10,:))';
    wav(:,:,j) = squeeze(q.wkmav(1,1:10,:))';
  end;
  
  % Plot rewards and reward predicaaamyfig(4,7.5);
  h = plot(rew(:,1),'linewidth',2); set(h,'color',colr1);
  hold on;
  h = plot(m,'linewidth',0.2); set(h,'color',colv1);
  set(gca,'ylim',[-2.4 2.4],'xlim',[0 180],'ytick',-2:2:2,'xtick',0:20:180,'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rewards_and_predictions.eps']);
  close(gcf);
  
  % Plot MBON firing rates
  myfig(2,7.5);
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
  
  % Plot synaptic weights
  myfig(4,7.5);
  for j=1:10
    plot(wap(:,:,j),'linewidth',0.2,'color',colm1);
    hold on;
    plot(wav(:,:,j),'linewidth',0.2,'color',colm2);
  end;
  set(gca,'ylim',[0 0.25],'xlim',[0 180],'ytick',0:0.1:0.2,'xtick',[0 20:20:180],'xticklabel',{0 [] [] 60 [] [] 120 [] [] 180},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'weights_ap_and_av.eps']);
  close(gcf);
  
  % Plot MBON firing rates for extended, noisy rewards, to show that they
  % are not unstable
  go = zeros(820,10);
  nogo = zeros(820,10);
  for k=1:10
    rews = cell(10,1);
    for j=1:10
      rews{j} = mb_reward_schedules(13,(k-1)*10+j,0.1,180,2);
    end;
    rew = rews{1};
    for j=2:10
      rew = [rew; rews{j}(21:end,:)];
    end;
    q = mb_mv_a(1,k,820,rew,0.5*2.5e-2,'choose1',true);
    go(:,k) = q.go(:,1); nogo(:,k) = q.nogo(:,1);
  end;
  mgo = mean(go,2);
  mnogo = mean(nogo,2);
  myfig(2,7.5);
  h = plot(go,'linewidth',0.2); set(h,'color',colm1);
  hold on;
  h = plot(nogo,'linewidth',0.2); set(h,'color',colm2);
  set(gca,'ylim',[0 2.2],'xlim',[0 820],'ytick',0:2,'xtick',[0 200:200:820],'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_firing_rates_stability.eps']);
  close(gcf);
end;

%%%
%%% Compare MV-Eq.7 and VS-lambda (limited learning) models when
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
  %%%% Generate associative conditioning experiment data for MV-Eq.7 model
  %%%%
  %%%% To test the Eq. 8 version of the MV model, replace 'eq7' with 'eq8'
  %%%% on lines 532 and 569 
  
  % # RNG seeds
  nrep = 1000;
%     % Matlab R2012 
%     % MacBook
%     % Setup parallel pool
%     p = parcluster;
%     NC = min(nrep,p.NumWorkers - 1);
%     p = matlabpool('size');
%     if p<1
%       matlabpool(NC);
%     else
%       matlabpool close force local;
%       matlabpool(NC);
%     end;

  % Matlab R2018
  % Setup parallel pool
  p = gcp('nocreate');
  NC = min(nrep,10);
  if isempty(p)
    pp = parpool(NC);
  else
    if p.NumWorkers~=NC
      delete(p);
      pp = parpool(NC);
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
  genetic_int = [1 0];
  
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
              q = mb_mv_conditioning((it-1)*NC+labindex,rew,1.25e-2,'intervene_id',intrvn_id(targ),genetic_int(gi),intrvn_schedule,'plasticity_rule','eq7');
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
        q = mb_mv_conditioning((it-1)*NC+labindex,rew,1.25e-2);
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
  save([fpath 'associative_conditioning_data_for_MV_model_eq7.mat'],'pre','dec','cpre','cdec');
end;

if flag==9
  %%%%
  %%%% Generate associative conditioning experiment data for VS-lambda model
  %%%%
  
  % # RNG seeds
  nrep = 1000;
  % Setup parallel pool
  %     % Matlab R2012
  %     % MacBook
  %     % Setup parallel pool
  %     p = parcluster;
  %     NC = min(nrep,p.NumWorkers - 1);
  %     p = matlabpool('size');
  %     if p<1
  %       matlabpool(NC);
  %     else
  %       matlabpool close force local;
  %       matlabpool(NC);
  %     end;
  % Matlab R2018
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
  genetic_int = [1 0];
  
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
              q = mb_vs_conditioning((it-1)*NC+labindex,rew,5e-2,'intervene_id',intrvn_id(targ),genetic_int(gi),intrvn_schedule);
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
  %%%% Fig. 5d-e, analyse VS_lambda and MV-Eq.8 model learned behaviours in
  %%%% different experimental interventions compared to controls in an associative
  %%%% conditioning paradigm. 
  %%%%
  
  m = cell(2,1);  
  m{1} = load([fpath 'associative_conditioning_data_for_MV_model_eq8.mat']);  
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
% 	z = xlsread([fpath 'supptable2.xlsx'],'A3:AG167');
  z = xlsread([fpath 'source_data_1.xlsx'],'PIs_all_experiments','A3:AH167');

  testtimes = z(:,end);
%   z = z(:,[1 32]);
  z = z(:,[1 31]);
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
  
  % List specific examples to be used later
  exampleind = [33 53 85 89, 1 8 22]; 
  exind=[];
  for k=1:length(exampleind)
    ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
    exind = [exind; find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi))];
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
        ind = ind(~ismember(ind,exind));
        for k=1:length(ind)
          for q=1:nbl
            plot(plmod{j}(q,ind(k)),plexp{j}(q,ind(k))+0.25*(rand-.5),'o','color',colneurons{m},'markersize',1+2*weights(q,ind(k)));
            hold on;
          end;
        end;
      end;
    end;    
    % Plot weighted linear fit
    x = [-6.6 6.6]; plot(x,b'*[1,1;-6.6,6.6],'linewidth',1,'color',[0.5 0.5 0.5]);
%     set(gca,'ylim',[-14 10],'xlim',[-6.6 6.6],'ytick',-10:5:10,'xtick',-6:3:6,'box','off');
    set(gca,'ylim',[-3.1 2.6],'xlim',[-7.1 7.8],'ytick',-2:1:2,'xtick',-6:3:6,'box','off');
    set(gca,'fontname','helvetica','fontsize',6);
    hold off; myticks;
    
    % Print out Pearson correlation coefficient
    pcorr = corr(plmod{j}(:).*s.w,plexp{j}(:).*s.w); 
    if j==1
      fprintf(['MV Model: Corr coef Robust Fit = %f\n'],pcorr);
    elseif j==2
      fprintf(['VS Model: Corr coef Robust Fit = %f\n'],pcorr);
    end;
    
    % Compute p-value for Pearson correlation using bootstrapping (10^4 resamples)
    rp = zeros(1e4,1); % corr-coeffs to compute p-values using permutation test
    rc = zeros(1e4,1); % corr-coeffs to compute confidence intervals using bootstrapping
    for k=1:1e4 
      % p-value
      rp1 = randperm(length(plmod{j}(:))); % Random permutation indeces into model data 
      [b s] = robustfit(plmod{j}(rp1),plexp{j}(:)); %  New robustfit on randomly permuted data
      rp(k) = corr(plmod{j}(rp1)'.*s.w,plexp{j}(:).*s.w); % Pearson correaltion on randomly permuted data
      % Confidence interval
      rp1 = randi(length(plmod{j}(:)),length(plmod{j}(:)),1);
      [b s] = robustfit(plmod{j}(rp1),plexp{j}(rp1)); %  New robustfit on randomly permuted data
      rc(k) = corr(plmod{j}(rp1).*s.w,plexp{j}(rp1).*s.w); % Pearson correaltion on randomly permuted data
    end;
    pval = sum(rp>pcorr)/1e4;
    cinterval = quantile(rc,[0.05 0.95]);
    if j==1
      fprintf(['MV Model: p-value = %f\n'],pval);
      fprintf(['MV Model: confidence interval = %f %f\n'],cinterval(1),cinterval(2));
    elseif j==2
      fprintf(['VS Model: p-value = %f\n'],pval);
      fprintf(['VS Model: confidence interval = %f %f\n'],cinterval(1),cinterval(2));
    end;
   
    if j==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_dap_dav_vs_exp_pi_diff.eps']);
    elseif j==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_decm_da_vs_exp_pi_diff.eps']);
    end;
    delete(gcf);
  end;
  
  % Plot examples from above that will be highlighted in the text
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
          plot(plmod{j}(q,ind(m)),plexp{j}(q,ind(m))+0.2*(rand-.5),'o','color',colneurons{targ},'markersize',1+2*weights(q,ind(m)));
          hold on;
        end;
      end;
    end;hold off;    
%     set(gca,'ylim',[-14 10],'xlim',[-6.6 6.6],'ytick',-10:5:10,'xtick',-6:3:6,'box','off');
    set(gca,'ylim',[-3.1 2.6],'xlim',[-7.1 7.8],'ytick',-2:1:2,'xtick',-6:3:6,'box','off');
    myticks;
    if j==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_dap_dav_vs_exp_examples.eps']);
    elseif j==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_decm_da_vs_exp_examples.eps']);
    end;
    % Plot second set of examples, coloured by neuron type.
    for k=[7 6 5]
      ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
      ind = find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi));
      for q=1:nbl
        for m=1:length(ind)
          plot(plmod{j}(q,ind(m)),plexp{j}(q,ind(m))+0.2*(rand-.5),'o','color',(colneurons{targ,1}-0.05).^2,'markersize',1+2*weights(q,ind(m)));
          hold on;
        end;
      end;
    end; hold off;
%     set(gca,'ylim',[-14 10],'xlim',[-6.6 6.6],'ytick',-10:5:10,'xtick',-6:3:6,'box','off');
    set(gca,'ylim',[-3.1 2.6],'xlim',[-7.1 7.8],'ytick',-2:1:2,'xtick',-6:3:6,'box','off');
    myticks;
    if j==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_dap_dav_vs_exp_examples2.eps']);
    elseif j==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_decm_da_vs_exp_examples2.eps']);
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
  %%%% Supplementary Fig. 6a-c, testing VS_lambda and MV model performance 
  %%%% as a function of the number of cues (from 2-200).
  %%%%
  
  %%% Setup parallel pool
  nrep = 100;
%   % MATLAB R2012
%   p = parcluster;
%   NC = min(nrep,p.NumWorkers - 1);
%   p = matlabpool('size');
%   if p<1
%     matlabpool(NC);
%   else
%     matlabpool close force local;
%     matlabpool(NC);
%   end;
  % MATLAB R2018
  % Setup parallel pool
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
  percmax = zeros(length(nstim),nrep); % Percentage of trials in which the maximally rewarding cue was chosen
  prewbeta5 = zeros(length(nstim),nrep); % Probabilities of obtaining the reward for 100% accurate agent, beta = 5.
  rmse = zeros(length(nstim),nrep); % Root mean squared error between predictions and reinforcements
  nit = floor(nrep/NC) + 1; % # iterations of spmd
  tic;
% % %   for model=1:4
  for model=3
    for no=1:length(nstim)
      for it=1:nit
        spmd(min(NC,nrep-(it-1)*NC))
          % Reward schedule
          rew = mb_reward_schedules(3,(it-1)*NC+labindex,0,250,nstim(no),2,10);
          rew(1:50,:) = [];
          
          prew = exp(rew/0.2) ./ repmat(sum(exp(rew / 0.2),2),1,nstim(no));
          mrew = sum(prew .* rew,2); % Mean reward (reward x probabiity of choosing rewarding cue (from softmax))
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
            q = mb_vs(1,(it-1)*NC+labindex,200,rew,1e-1,'no',nstim(no),'lambda',12,'nk',nstim(no)*100);
          end;
        end;
        for j=1:min(NC,nrep-(it-1)*NC)
          qq = q{j};
          if model==1
            rewpercue1(no,(it-1)*NC+j) = sum(max(rew{j},[],2)) / 200; % Max reward per trial
            rewpercue2(no,(it-1)*NC+j) = sum(mrew{j}) / 200; % Weighted mean reward per trial
            rewpercue3(no,(it-1)*NC+j) = sum(mean(qq.r,2)) / 200; % Mean reward per trial
%             prewbeta5(no,(it-1)*NC+j) = mean(max(prew{j},[],2));
          end;
          if model==3
            [dum mxind] = find(rew{j}==(max(rew{j},[],2)*ones(1,nstim(no))));
            [s sind] = sort(dum); clear dum s;
            mxind = mxind(sind);
            percmax(no,(it-1)*NC+j) = mean(mxind == qq.decision) * 100;
            rmse(no,(it-1)*NC+j) = sqrt(mean(mean((rew{j} - (qq.go - qq.nogo)).^2,2),1));
            prewbeta5(no,(it-1)*NC+j) = mean(max(prew{j},[],2));
          end;          
          obr = qq.sr; % Obtained reward per trial          
          rewpercue4{model}(no,(it-1)*NC+j) = obr / 200;
        end;
        clear rew q mrew prew;
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
  myfig(2,7.5);
  h = fill([nstim'; fliplr(nstim)'],[mean(prewbeta5,2)+std(prewbeta5,1,2); flipud(mean(prewbeta5,2)-std(prewbeta5,1,2))],'k');
  set(h,'edgecolor','none','facealpha',0.35); hold on;
  h = plot(nstim',mean(prewbeta5,2),'linewidth',1,'color','k');
  h = fill([nstim'; fliplr(nstim)'],[mean(percmax/100,2)+std(percmax/100,1,2); flipud(mean(percmax/100,2)-std(percmax/100,1,2))],'k');
  set(h,'edgecolor','none','facealpha',0.35); 
  h = plot(nstim',mean(percmax/100,2),'--','linewidth',1,'color','k');
  set(gca,'xscale','log','ylim',[0 1],'ytick',0:0.5:1.5,'xlim',[2 200],'xtick',[2:2:10 20:20:100 200],'xticklabel',{2 [] [] [] 10 20 [] [] [] 100 200},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'prob_max_reward.eps']); 
  
  h = fill([nstim'; fliplr(nstim)'],[mean(rmse,2)+std(rmse,1,2); flipud(mean(rmse,2)-std(rmse,1,2))],'k');
  set(h,'edgecolor','none','facealpha',0.35); hold on;
  h = plot(nstim',mean(rmse,2),'linewidth',1,'color','k');
  set(gca,'xscale','log','ylim',[0.05 0.45],'ytick',[0.1 0.4],'xlim',[2 200],'xtick',[2:2:10 20:20:100 200],'xticklabel',{2 [] [] [] 10 20 [] [] [] 100 200},'box','off');
  myticks; hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'rmse.eps']);
  delete(gcf);
end;

if flag==14
  %%%%
  %%%% Supplementary Fig. 6d, testing differences in performance between
  %%%% the perfect plasticity agent and the MV model as a function of
  %%%% noise in the rewards and the # of cues.
  %%%%
  
  %%% Setup parallel pool
  nrep = 100;
  % Matlab R2012
%   p = parcluster;
%   NC = min(nrep,p.NumWorkers - 1);
%   p = matlabpool('size');
%   if p<1
%     matlabpool(NC);
%   else
%     matlabpool close force local;
%     matlabpool(NC);
%   end;
      % MATLAB R2018
      % Setup parallel pool
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
  %%%% Supplementary Fig. 6e-g, Comparing reward predictions for single
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

if flag==16
  %%%%
  %%%% Supp Fig. 7, comparing different models. 
  %%%% 1) Compare VSlambda against MV-Eq.7 model
  %%%% 2) Compare MV-Eq.7 model (e.g. dw+ ~ (d+ - d-)) against MV-Eq.8 model 
  %%%%    (e.g. dw+ ~ (w_k*k - d-)) 
  %%%% 3) Compare the MV-Eq.8 model agains the experimental data.
  %%%%
  
  m = cell(3,1);
  m{1} = load([fpath 'associative_conditioning_data_for_VS_model.mat']);
  m{2} = load([fpath 'associative_conditioning_data_for_MV_model_eq7.mat']);
  m{3} = load([fpath 'associative_conditioning_data_for_MV_model_eq8.mat']);
  nrep = 1000; % # repeat runs of the simulation per condition
  nip = 4; % # intervention protocols (intervene during training/testing/both)
  ntarg = 4; % # target neurons
  ngi = 2; % # genetic interventions (Shibire/TrpA)
  nrt = 3; % # reward types (+ve/-ve)
  nt = 30; % # trials
  nav = 50; % # runs over each block of averages
  nbl = nrep / nav; % # blocks of averages
  nmod = 3; % # models
  cc = cell(3,1); % Cell array containg performance index (PI) data from simulated genetic intervention conditions
	ct = cell(3,1); % Cell array containg performance index (PI) data from simulated control conditions
  plmod = cell(3,1); % Cell array for plotting data from the model
  plexp = cell(1,1); % Cell array for plotting data from experiments
  
  %%% Orgainse model PI data (experimental conditions) for further analysis
  for model=1:3
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
	for model=1:3
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
  z = xlsread([fpath 'source_data_1.xlsx'],'PIs_all_experiments','A3:AH167');

  testtimes = z(:,end);
  z = z(:,[1 31]);
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
  
  % Compute test statistics for all simulations in all models
  ts = cell(3,1);
  for model=1:3
    ts{model} = zeros(nbl,size(condslist,1)*nbl);
    for j=1:size(condslist,1)
      if condslist(j,4)==1
        temp = ct{model}(1:2:end,1);
        cont = mean(temp(:))*ones(nbl,1);
      elseif condslist(j,4)==2
        temp = ct{model}(2:2:end,1);
        cont = mean(temp(:))*ones(nbl,1);
      elseif condslist(j,4)==3
        cont = zeros(nbl,1);ones(nbl,1);
      end;
      
      temp = (cc{model}(indcond==j,1) + 1) / 2;
      cont = (cont + 1) / 2;
      ts{model}(:,j) = (temp - cont) ./ sqrt((temp + cont)/2.*(1 - (temp+cont)/2)*2/nav);
    end;
  end;
  
  % Populate the cell arrays with test statistics for plotting against
  % experimental test statistics
  for model=1:3
    if model==1
      plexp{model} = zeros(nbl,length(ind0));
    end;
    plmod{model} = zeros(nbl,length(ind0));
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
      if model==1
        % Make repeat copies of each experimetnal data point, to pair up
        % with the 20 model data points for each condition
        plexp{model}(:,j) = z(ind0(j),2) * ones(nbl,1);
      end;
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
  
  %%%
  %%% 1) Plot VSlambda versus MV-Eq.8 test statistics
  %%%
  myfig(8,8);
  % Plot, coloured by neuron type
  for j=1:size(condslist,1)
    plot(ts{3}(:,j)+0.2*(rand(nbl,1) - 0.5),ts{1}(:,j)+0.2*(rand(nbl,1) - 0.5),'o','color',colneurons{condslist(j,2)},'markersize',2);
    hold on;
    if (abs(mean(ts{3}(:,j)))<1.5 && abs(mean(ts{1}(:,j)))>3) || (abs(mean(ts{1}(:,j)))<1.5 && abs(mean(ts{3}(:,j)))>5)
      fprintf('cond=%d%d%d%d\n',condslist(j,:));
    end;
  end;
  set(gca,'ylim',[-10 10],'xlim',[-10 10],'ytick',-10:5:10,'xtick',-10:5:10,'box','off');
  hold off;
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_VS_vs_MV8_pi_diff.eps']);
  delete(gcf);
  
  %%%
  %%% 2) Plot MV-Eq.7 versus MV-Eq.8 test statistics
  %%%
  % Examples in condslist to point out in figure, which correspond to the
  % conditions:
  % condslist(18,:) = [1 3 2 3] -- includes experimental evidence in supptable2
  % condslist(24,:) = [1 4 2 3] -- includes experimental evidence in supptable2
  % condslist(88,:) = [4 3 2 1]
  % condslist(95,:) = [4 4 2 2]
  exampleind = [18, 24, 88, 95];
  myfig(8,8);
  % Plot, coloured by neuron type
  ind = 1:size(condslist,1); ind(exampleind) = [];
  for j=ind
    plot(ts{3}(:,j)+0.2*(rand(nbl,1) - 0.5),ts{2}(:,j)+0.2*(rand(nbl,1) - 0.5),'o','color',colneurons{condslist(j,2)},'markersize',2);
    hold on;
    if (abs(mean(ts{3}(:,j)))<1.5 && abs(mean(ts{2}(:,j)))>3) || (abs(mean(ts{2}(:,j)))<1.5 && abs(mean(ts{3}(:,j)))>5)
      fprintf('cond=%d%d%d%d\n',condslist(j,:));
    end;
  end;
  set(gca,'ylim',[-10 10],'xlim',[-10 10],'ytick',-10:5:10,'xtick',-10:5:10,'box','off');
  hold off;
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_MV8_vs_MV7_pi_diff.eps']);
  for j=exampleind
    plot(ts{3}(:,j)+0.2*(rand(nbl,1) - 0.5),ts{2}(:,j)+0.2*(rand(nbl,1) - 0.5),'o','color',colneurons{condslist(j,2)},'markersize',2);
    hold on;
  end;
  set(gca,'ylim',[-10 10],'xlim',[-10 10],'ytick',-10:5:10,'xtick',-10:5:10,'box','off');
  hold off;
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_MV8_vs_MV7_pi_diff_examples.eps']);
  delete(gcf);
  
  %%%
  %%% 3) Plot the experiment versus MV-Eq.7 model test statistics
  %%%
  % List specific examples to be used later
  exampleind = [33 53 85 89, 8 22 1]; 
  exind=[];
  for k=1:length(exampleind)
    ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
    exind = [exind; find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi))];
  end;  
  myfig(8,8);
  % Compute the weighted linear fit
  [b s] = robustfit(plmod{2}(:),plexp{1}(:));
  % Weights from robustfit
  weights = reshape(s.w,nbl,length(ind0));
  % Plot, coloured by neuron type
  for m=[1 3 4 2]
    for n=1:2
      ind = find(ind0targ(:,m)&ind0gi(:,n));
      ind = ind(~ismember(ind,exind));
      for k=1:length(ind)
        for q=1:nbl
          plot(plmod{2}(q,ind(k)),plexp{1}(q,ind(k))+0.5*(rand-.5),'o','color',colneurons{m},'markersize',1+2*weights(q,ind(k)));
          hold on;
        end;
      end;
    end;
  end;
  % Plot weighted linear fit
  x = [-6.6 6.6]; plot(x,b'*[1,1;-6.6,6.6],'linewidth',1,'color',[0.5 0.5 0.5]);
  set(gca,'ylim',[-3.1 2.6],'xlim',[-7.1 7.8],'ytick',-2:2:2,'xtick',-6:6:6,'box','off');
  hold off; myticks;
  
  % Print out Pearson correlation coefficient
  pcorr = corr(plmod{2}(:).*s.w,plexp{1}(:).*s.w);
  fprintf(['MV7 Model: Corr coef Robust Fit = %f\n'],pcorr);
  
  % Compute p-value for Pearson correlation using bootstrapping (10^4 resamples)
  rp = zeros(1e4,1); % corr-coeffs to compute p-values using permutation test
  rc = zeros(1e4,1); % corr-coeffs to compute confidence intervals using bootstrapping
  for k=1:1e4
    % p-value
    rp1 = randperm(length(plmod{2}(:))); % Random permutation indeces into model data
    [b s] = robustfit(plmod{2}(rp1),plexp{1}(:)); %  New robustfit on randomly permuted data
    rp(k) = corr(plmod{2}(rp1)'.*s.w,plexp{1}(:).*s.w); % Pearson correaltion on randomly permuted data
    % Confidence interval
    rp1 = randi(length(plmod{2}(:)),length(plmod{2}(:)),1);
    [b s] = robustfit(plmod{2}(rp1),plexp{1}(rp1)); %  New robustfit on randomly permuted data
    rc(k) = corr(plmod{2}(rp1).*s.w,plexp{1}(rp1).*s.w); % Pearson correaltion on randomly permuted data
  end;
  pval = sum(rp>pcorr)/1e4;
  cinterval = quantile(rc,[0.05 0.95]);
  fprintf(['MV7 Model: p-value = %f\n'],pval);
  fprintf(['MV7 Model: confidence interval = %f %f\n'],cinterval(1),cinterval(2));

  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_MV7_vs_exp_pi_diff.eps']);
  delete(gcf);
  
  % Plot examples from above that will be highlighted in the text
  myfig(8,8);
  % Plot first set of examples, coloured by neuron type.
  for k=1:4
    ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
    ind = find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi));
    for m=1:length(ind)
      for q=1:nbl
        plot(plmod{3}(q,ind(m)),plexp{1}(q,ind(m))+0.2*(rand-.5),'o','color',colneurons{targ},'markersize',1+2*weights(q,ind(m)));
        hold on;
      end;
    end;
  end;hold off;
  set(gca,'ylim',[-3.1 2.6],'xlim',[-7.1 7.8],'ytick',-2:1:2,'xtick',-6:3:6,'box','off');
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_MV7_vs_exp_examples.eps']);
  % Plot second set of examples, coloured by neuron type.
  for k=5:7
    ip = ex(exampleind(k),1); targ = ex(exampleind(k),2); gi = ex(exampleind(k),3); rt = ex(exampleind(k),4);
    ind = find(ind0ip(:,ip)&ind0rt(:,rt)&ind0targ(:,targ)&ind0gi(:,gi));
    for q=1:nbl
      for m=1:length(ind)
        plot(plmod{3}(q,ind(m)),plexp{1}(q,ind(m))+0.2*(rand-.5),'o','color',(colneurons{targ}-0.05).^2,'markersize',1+2*weights(q,ind(m)));
        hold on;
      end;
    end;
  end; hold off;
  set(gca,'ylim',[-3.1 2.6],'xlim',[-7.1 7.8],'ytick',-2:1:2,'xtick',-6:3:6,'box','off');
  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_model_MV7_vs_exp_examples2.eps']);
  delete(gcf);
end;

if flag==17
  %%%%
  %%%% Supplementary Fig. 4, comparing the two MV model learning rules
  %%%%
  seed = 3;
  rew = zeros(30,2) + 0.1*randn(30,2); 
  go = zeros(10,10,2);
  nogo = zeros(10,10,2);
  mhat = zeros(10,10,2);
  for j=1:10
    % For Eq. 7 learning rule
    q = mb_mv_conditioning(j,rew,1.25e-2,'intervene_id',3,0,[5*ones(10,1);zeros(20,1)],'plasticity_rule','eq7');
    go(:,j,1) = q.go(1:10,1);
    nogo(:,j,1) = q.nogo(1:10,1);
    mhat(:,j,1) = q.go(1:10,1) - q.nogo(1:10,1);
    % For Eq. 8 learning rule
    q = mb_mv_conditioning(j,rew,2.5e-2,'intervene_id',3,0,[5*ones(10,1);zeros(20,1)],'plasticity_rule','eq8');
    go(:,j,2) = q.go(1:10,1);
    nogo(:,j,2) = q.nogo(1:10,1);
    mhat(:,j,2) = q.go(1:10,1) - q.nogo(1:10,1);
  end;
  
  % Plot reward predictions
  myfig(4,7.5);
  h = plot(mhat(:,:,1),'linewidth',0.5); set(h,'color',colv1);
  hold on;
  h = plot(mhat(:,:,2),'linewidth',0.5); set(h,'color',colv1.^(1/4));
  plot([1 10],[0 0],'k','linewidth',1);
  set(gca,'ylim',[-0.4 2.6],'xlim',[1 10],'ytick',0:2,'xtick',[1 10],'box','off');
  set(gca,'fontname','helvetica','fontsize',6); hold off; myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'predictions_both_plasticityrules.eps']);
  
  % Plot MBON firing rates for the CS+
  h = plot(go(:,:,1),'linewidth',0.5); set(h,'color',colm1);
  hold on;
  h = plot(nogo(:,:,1),'linewidth',0.5); set(h,'color',colm2);
  set(gca,'ylim',[0 2.6],'xlim',[1 10],'ytick',0:2,'xtick',[1 10],'box','off');
  set(gca,'fontname','helvetica','fontsize',6); hold off; myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_rates_eq7.eps']);
  h = plot(go(:,:,2),'linewidth',0.5); set(h,'color',colm1.^(1/4));
  hold on;
  h = plot(nogo(:,:,2),'linewidth',0.5); set(h,'color',colm2.^(1/4));
  set(gca,'ylim',[0 2.6],'xlim',[1 10],'ytick',0:2,'xtick',[1 10],'box','off');
  set(gca,'fontname','helvetica','fontsize',6); hold off; myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) 'mbon_rates_eq8.eps']);
  close(gcf);
end;

if flag==18
  %%%%
  %%%% Supplementary Fig. 8, comparing the VSlambda and MV models with the 
  %%%% experiments of Felsenberg et al. (2017 & 2018)
  %%%%
  
  % Model parameters
  nt = 40; % # time steps per simulation
  nrep = 1000; % # repeats per condition
  NC = min(nrep,10); % # threads for parallelisation
  nit = floor(nrep/NC) + 1; % # iterations of spmd
  nav = 50; % # simulations over which to average
  nbl = nrep / nav; % # blocks of averages (i.e. test stats per simulated condition)
  
  %%% Load model data if it already exists
  d = dir(fpath);
  for j=1:length(d)
    if strcmp(d(j).name,'felsenberg_extinction_data_for_both_models.mat')
      load([fpath 'felsenberg_extinction_data_for_both_models.mat']);
    end;
  end;
  
  %%% Setup parallel pool
  if ~exist('dec','var')
    p = gcp('nocreate');
    if isempty(p)
      pp = parpool(NC);
    else
      if p.NumWorkers~=NC
        delete(p);
        pp = parpool(NC);
      end;
    end;
  end;
  
  % Experiment parameters
  reexp = [1 2]; % Stimulus ID for re-exposure
  targ = [1 2 3 4]; % Target cell ID for blocking: 0 - none; 1 - M+; 2 - M-; 3 - D+; 4 - D-
  rtype = [-1 1]; % Reinforcement type: punishment (-1) or reward (+1)
  sched = [ones(20,1); 0.1*ones(10,1); ones(10,1)]; % Intervention schedule: block during re-exposure only
  
  % Condition codes
  % [A B C]: A - Target neuron: M+ (1), M- (2), D+ (3), D- (4) 
  %          B - CS+ (1), CS- (2)
  %          C - Punishment (1), Reward (2)
  conds = [1 1 2;...
           1 2 2;...
           2 1 1;...
           2 2 1;...
           3 1 1;...
           3 1 2;...
           3 2 1;...
           3 2 2;...
           4 1 1;...
           4 1 2;...
           4 2 1;...
           4 2 2];
  
  %%% Run simulations
  if ~exist('dec','var')
    % Generate data for conditions
    dec = cell(2,1); % For storing test decisions
    for j=1:2
      dec{j} = zeros(size(conds,1),nrep,2);
    end;
    for ii=1:size(conds,1)
      iTarg = conds(ii,1);
      iReexp = conds(ii,2);
      iRtype = conds(ii,3);
      for it=1:nit
        spmd(min(NC,nrep-(it-1)*NC))
          % Contruct reward schedule
          rew = zeros(nt,2);
          rew(1:10,1) = rtype(iRtype);
          rew = rew + 0.1*randn(nt,2);
          % Run model
          qvs = mb_vs_conditioning_fels((it-1)*NC+labindex,rew,5e-2,'intervene_id',iTarg,1,sched,'reexp_id',iReexp);
          qmv = mb_mv_conditioning_fels((it-1)*NC+labindex,rew,1.25e-2,'intervene_id',iTarg,1,sched,'reexp_id',iReexp);
        end;
        for j=1:min(NC,nrep-(it-1)*NC)
          qqvs = qvs{j};
          qqmv = qmv{j};
          
          dec{1}(ii,(it-1)*NC+j,:) = qqvs.decision(31:32);
          dec{2}(ii,(it-1)*NC+j,:) = qqmv.decision(31:32);
        end;
        clear rew qvs qmv;
      end;
    end;
    
    % Generate data for controls
    cdec = cell(2,1); % For storing test decisions
    for j=1:2
      cdec{j} = zeros(2,2,nrep,2);
    end;
    for iReexp=1:2
      for iRtype=1:2
        for it=1:nit
          spmd(min(NC,nrep-(it-1)*NC))
            % Contruct reward schedule
            rew = zeros(nt,2);
            rew(1:10,1) = rtype(iRtype);
            rew = rew + 0.1*randn(nt,2);
            % Run model
            qvs = mb_vs_conditioning_fels((it-1)*NC+labindex,rew,5e-2,'reexp_id',iReexp);
            qmv = mb_mv_conditioning_fels((it-1)*NC+labindex,rew,1.25e-2,'reexp_id',iReexp);
          end;
          for j=1:min(NC,nrep-(it-1)*NC)
            qqvs = qvs{j};
            qqmv = qmv{j};
            
            cdec{1}(iReexp,iRtype,(it-1)*NC+j,:) = qqvs.decision(31:32);
            cdec{2}(iReexp,iRtype,(it-1)*NC+j,:) = qqmv.decision(31:32);
          end;
          clear rew qvs qmv;
        end;
      end;
    end;
    save([fpath 'felsenberg_extinction_data_for_both_models.mat'],'dec','cdec');
  end;
  
  %%% Compute test statistics from simulation data
  % Init arrays for test statistic
  plmod = cell(2,1);
  for j=1:2
    plmod{j} = zeros(size(conds,1),nbl);
  end;
  
  cont = cell(2,1); % Control simulations for each model
  for j=1:2
    cont{j} = zeros(2,2,nbl);
  end;
  
  % Compute fractions of CS+ choices for each condition and control
  for model=1:2
    for iReexp=1:2
      for iRtype=1:2
        for ibl=1:nbl
          n1 = sum(sum(cdec{model}(iReexp,iRtype,(ibl-1)*nav+1:ibl*nav,:)==1,3),4); % # CS+ choices
          cont{model}(iReexp,iRtype,ibl) = n1 / (numel(cdec{model}(iReexp,iRtype,(ibl-1)*nav+1:ibl*nav,:))); % Fraction of CS+ choices
        end; 
      end;
    end;
    
    for ii=1:size(conds,1)
      iReexp = conds(ii,2);
      iRtype = conds(ii,3);
      for ibl=1:nbl
        n1 = sum(sum(dec{model}(ii,(ibl-1)*nav+1:ibl*nav,:)==1,2),3); % # CS+ choices
        temp = n1 / (numel(dec{model}(ii,(ibl-1)*nav+1:ibl*nav,:))); % Fraction of CS+ choices
        plmod{model}(ii,ibl) = (temp - cont{model}(iReexp,iRtype,ibl)) ./ sqrt((temp + cont{model}(iReexp,iRtype,ibl))/2.*(1 - (temp+cont{model}(iReexp,iRtype,ibl))/2)*2/nav);        
      end;
      indnan = isnan(plmod{model}(ii,:));
      plmod{model}(ii,indnan) = 0;
    end;
  end;
 
  %%% Import Felsenberg et al. data
  z = xlsread([fpath 'source_data_2.xlsx'],'PIs_Felsenberg_et_al','A3:AF14');
  plexp = z(:,15);
  exp4corr = reshape(plexp*ones(1,nbl),length(plexp)*nbl,1); % For computing correlations
  
  %%% Plot experimental versus model data
  % Generate cell array of colours for plotting
  colneurons = [colm1; colm2; coldap; coldav];
  
  % Plot data
  for model=1:2
    % Compute the weighted linear fit
    [b s] = robustfit(plmod{model}(:),exp4corr);       
    % Plot, coloured by neuron type   
    myfig(8,8);
    for ii=1:size(conds,1)
      iTarg = conds(ii,1);
      plot(squeeze(plmod{model}(ii,:)),plexp(ii)+0.25*(rand(nbl,1)-.5),'o','color',colneurons(iTarg,:),'markersize',3);
      hold on;
    end;
    set(gca,'!(ylim',[-1.3 1.7],'xlim',[-10.1 10.1],'ytick',-1:0.5:1.5,'xtick',-10:5:10,'box','off');
    set(gca,'fontname','helvetica','fontsize',6); hold off; myticks;
    
    % Print out Pearson correlation coefficient
    pcorr = corr(plmod{model}(:).*s.w,exp4corr.*s.w);
    if model==1
      fprintf(['VS Model: Corr coef Robust Fit = %f\n'],pcorr);
    elseif model==2
      fprintf(['MV8 Model: Corr coef Robust Fit = %f\n'],pcorr);
    end;
    
    % Compute p-value for Pearson correlation using bootstrapping (10^4 resamples)
    rp = zeros(1e4,1); % corr-coeffs to compute p-values using permutation test
    rc = zeros(1e4,1); % corr-coeffs to compute confidence intervals using bootstrapping
    for k=1:1e4
      % p-value
      rp1 = randperm(length(plmod{model}(:))); % Random permutation indeces into model data
      [b s] = robustfit(plmod{model}(rp1),exp4corr); %  New robustfit on randomly permuted data
      rp(k) = corr(plmod{model}(rp1)'.*s.w,exp4corr.*s.w); % Pearson correaltion on randomly permuted data
      % Confidence interval
      rp1 = randi(length(plmod{model}(:)),length(plmod{model}(:)),1);
      [b s] = robustfit(plmod{model}(rp1),exp4corr(rp1)); %  New robustfit on randomly permuted data
      rc(k) = corr(plmod{model}(rp1).*s.w,exp4corr(rp1).*s.w); % Pearson correaltion on randomly permuted data
    end;
    pval = sum(rp>pcorr)/1e4;
    cinterval = quantile(rc,[0.05 0.95]);
    if model==1
      fprintf(['VS Model: p-value = %f\n'],pval);
      fprintf(['VS Model: confidence interval = %f %f\n'],cinterval(1),cinterval(2));
    elseif model==2
      fprintf(['MV8 Model: p-value = %f\n'],pval);
      fprintf(['MV8 Model: confidence interval = %f %f\n'],cinterval(1),cinterval(2));
    end;
    
    if model==1
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_VSmodel_vs_exp_pi_diff.eps']);
    elseif model==2
      myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_rho_MVmodel_vs_exp_pi_diff.eps']);
    end;
    delete(gcf);
  end;
  
elseif flag==19
  %%%
  %%% Fig. 6, blocking experiments for independent (variables with 1) and dependent
  %%% (variables with 0) stimulus representations.
  %%%
  m1 = zeros(30,2,50); m0 = zeros(30,2,50);
  d1 = zeros(50,20); d0 = zeros(50,20);
  for i1=1:50
    for i2=1:20
      q = mb_mv_blocking((i1-1)*20+i2,1.25e-2,'pflip',[0 0]);
      if i2==1
        m1(:,1,i1) = q.go(:,1)-q.nogo(:,1);
        m1(:,2,i1) = q.go(:,2)-q.nogo(:,2);
      end;
      d1(i1,i2) = mean(q.decision(21:22));
      q = mb_mv_blocking((i1-1)*20+i2,1.25e-2,'pflip',[0.8 0.2]);
      if i2==1
        m0(:,1,i1) = q.go(:,1)-q.nogo(:,1);
        m0(:,2,i1) = q.go(:,2)-q.nogo(:,2);
      end;
      d0(i1,i2) = mean(q.decision(21:22));
    end;
  end;
  
  mm1 = mean(mean(m1(21:22,:,:),1),3);
  mm0 = mean(mean(m0(21:22,:,:),1),3);
  sm1 = std(mean(m1(21:22,:,:),1),[],3);
  sm0 = std(mean(m0(21:22,:,:),1),[],3);
  
  d1 = mean(d1,1) - 1;
  d0 = mean(d0,1) - 1;
  d1 = 2 * d1 - 1;
  d0 = 2 * d0 - 1;
  md1 = mean(d1);
  md0 = mean(d0);
  sd1 = std(d1);
  sd0 = std(d0);
  
  myfig(5,2);
  % Independent reward prediction
  am = squeeze(mean(m1(21:22,:,:),1));
  g1 = exp(-(am(1,:) - mm1(1))'.^2/2/sm1(1)^2);
  g2 = exp(-(am(2,:) - mm1(2))'.^2/2/sm1(2)^2);
  mybar([1 2],[mm1(1) mm1(2)],'style','categorical','errorbar',sm1,'color',[0.5 0.5 0.5]); hold on;
  plot(ones(50,1)+(2*(rand(50,1)>0.5)-1).*(0.35*g1.*rand(50,1)),am(1,:)','ok','markersize',2);
  plot(2*ones(50,1)+(2*(rand(50,1)>0.5)-1).*(0.35*g2.*rand(50,1)),am(2,:)','ok','markersize',2);
  set(gca,'ylim',[-0.3 1.3],'xlim',[0.5 2.5],'ytick',-0.4:0.4:1.2,'xtick',[1 2],'box','off');  myticks;
  hold off;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_reward_predictions_indy1.eps']);
 
  % Non-independent reward prediction
  am = squeeze(mean(m0(21:22,:,:),1));
  g1 = exp(-(am(1,:) - mm0(1))'.^2/2/sm0(1)^2);
  g2 = exp(-(am(2,:) - mm0(2))'.^2/2/sm0(2)^2);
  mybar([1 2],[mm0(1) mm0(2)],'style','categorical','errorbar',sm0,'color',[0.5 0.5 0.5]); hold on;
  plot(ones(50,1)+(2*(rand(50,1)>0.5)-1).*(0.35*g1.*rand(50,1)),am(1,:)','ok','markersize',2);
  plot(2*ones(50,1)+(2*(rand(50,1)>0.5)-1).*(0.35*g2.*rand(50,1)),am(2,:)','ok','markersize',2);
  hold off;
  set(gca,'ylim',[-0.3 1.3],'xlim',[0.5 2.5],'ytick',-0.4:0.4:1.2,'xtick',[1 2],'box','off');  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_reward_predictions_indy0.eps']);
  
  % Performance index
  g1 = exp(-(d1 - md1).^2/2/sd1^2);
  g2 = exp(-(d0 - md0).^2/2/sd0^2);
  mybar([1 2],[md1 md0],'style','categorical','errorbar',[sd1 sd0],'color',[0.5 0.5 0.5]); hold on;
  plot(ones(1,20)+(2*(rand(1,20)>0.5)-1).*(0.35*g1.*rand(1,20)),d1,'ok','markersize',2);
  plot(2*ones(1,20)+(2*(rand(1,20)>0.5)-1).*(0.35*g2.*rand(1,20)),d0,'ok','markersize',2);
  hold off;
  set(gca,'ylim',[-0.15 0.7],'xlim',[0.5 2.5],'ytick',0:0.3:0.6,'xtick',[1 2],'box','off');  myticks;
  myprint(gcf,'epsc',[fpath 'flag' num2str(flag) '_performance_indeces_indy10.eps']);
end;


