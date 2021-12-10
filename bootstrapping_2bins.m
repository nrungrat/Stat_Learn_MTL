% This script computes resampling stats across conditions
%-----------------------------------------------------------------
close all; clear all;
tic
%load('controls_data_2bins_manuscriptfig.mat')
%load('patients_data_2bins_manuscriptfig.mat')

%load('controls_data_2bins_022519.mat')
%load('patients_data_2bins_022519.mat')

%load('controls_data.mat')
load('patients_data.mat')

numboot = 1000;

%color stuff
c_foc = [0 0.4 0]; %green
c_div = [0.4 0 0.6]; %purple
c_exp = [0 0.4 0.8];
c_neu = [0.5 0.5 0.5]; %grey
c_unexp = [0.8 0 0]; %orange
c_hi = [0.6 0 0.2]; %maroon
c_lo = [0.4 0.4 0.4]; %grey
diffcolor1 = [0.5 0.5 0.5]; %exp - unexp
diffcolor2 = [0.5 0.5 0.5];
opac = 0.6; %opacity value for plotting CI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Resp Trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1A) resp trajectories as a function of expectation
% bin 1
for nboot = 1:numboot
    %nboot
    fprintf('resp traj vs exp_bin1: %s\n', num2str(nboot))
    sub_ind = 1:numel(allsub);
    %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
    for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
        clear nexp_cond nun_cond ndiff_expun_avg %doing this b/c #of trials in each bin vary across subj
        
        exp_cond = data.dist_exp1{sub_ind(ss)}; %all exp trials in b1 of that subj
        un_cond = data.dist_un1{sub_ind(ss)}; %all un trials in b1 of that subj
        
        %%% 2 ways of doing this i think: 1) fixed # across subj determined by the
        %%% min # of un trials across subj; 2) the number of un of each subj
        %%% trying 2) first
        ind = randi(size(un_cond, 2), [1, size(un_cond, 2)]); % this is an index matrix to sample
        % out XX trials (with replacement) each from expected and
        % unexpected condition
        for timepoint = 1:frm
            sh_exp = randsample(size(exp_cond, 2), size(exp_cond, 2))'; %shuffle trials for bootstrapping
            exp_ind = sh_exp(ind);
            nexp_cond(timepoint, :) = exp_cond(timepoint, exp_ind(1, :)); %resampled expected trials
            nexp_avg(timepoint, :) = nanmean(nexp_cond(timepoint, :), 2);
            nun_cond(timepoint, :) = un_cond(timepoint, ind(1, :)); %resampled unexpected trials
            nun_avg(timepoint, :) = nanmean(nun_cond(timepoint, :), 2);
            ndiff_expun_avg(timepoint, :) = nun_avg(timepoint, :) - nexp_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
        end
        nexp_cond_b1{ss} = nexp_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        nun_cond_b1{ss} = nun_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nexp_avg_b1(:, ss) = nexp_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        nun_avg_b1(:, ss) = nun_avg;
        ndiff_expun_avg_b1(:, ss) = ndiff_expun_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    traj.nexp_savg_b1(:, nboot) = nanmean(nexp_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
    traj.nun_savg_b1(:, nboot) = nanmean(nun_avg_b1, 2);
    traj.ndiff_expun_savg_b1(:, nboot) = nanmean(ndiff_expun_avg_b1, 2);
end

for tp = 1:frm %loop over timepoints to find average across iterations
    %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    traj.ndiff_expun_b1_CI(:, tp) = prctile(traj.ndiff_expun_savg_b1(tp, :), [2.5 97.5]);
    traj.nexp_b1_CI(:, tp) = prctile(traj.nexp_savg_b1(tp, :), [2.5 97.5]);
    traj.nun_b1_CI(:, tp) = prctile(traj.nun_savg_b1(tp, :), [2.5 97.5]);
    traj.p_exp_b1(:, tp) = min((2*min(numel(find(traj.ndiff_expun_savg_b1(tp, :) > 0))./numboot, ...
        1-(numel(find(traj.ndiff_expun_savg_b1(tp, :) < 0))./numboot))), ...
        (2*min(numel(find(traj.ndiff_expun_savg_b1(tp, :) < 0))./numboot, ...
        1-(numel(find(traj.ndiff_expun_savg_b1(tp, :) > 0))./numboot))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1B) plotting traj as a function of exp (bin1-bin4)
for fig = 1
    figure(1); clf;
    n1 = suptitle(['Resp Traj as a f(n) of Expectation (n=' num2str(size(allsub,2)) ')'])
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    %bin1
    subplot(3, 1, 1) %resampling data of exp vs unexp
    p1 = plot(timex, nanmean(traj.nexp_savg_b1, 2), 'Color', c_exp, 'LineWidth', 3); hold on; %resampled expected trials
    p2 = plot(timex, traj.nexp_b1_CI(1, :), 'Color', c_exp, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, traj.nexp_b1_CI(2, :), 'Color', c_exp, 'LineWidth', 1.5); hold on;
    p2.Color(4) = opac;
    p3.Color(4) = opac;
    
    p4 = plot(timex, nanmean(traj.nun_savg_b1, 2), 'Color', c_unexp, 'LineWidth', 3); %resampled unexpected trials
    p5 = plot(timex, traj.nun_b1_CI(1, :), 'Color', c_unexp, 'LineWidth', 1.5); hold on;
    p6 = plot(timex, traj.nun_b1_CI(2, :), 'Color', c_unexp, 'LineWidth', 1.5);
    p5.Color(4) = opac;
    p6.Color(4) = opac;
    
    ylabel('Distance (a.u.)')
    legend([p1 p4], {'Expected', 'Unexpected'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %unexpected minus expected
    plot(timex, nanmean(traj.ndiff_expun_savg_b1, 2), 'Color', diffcolor1, 'LineWidth', 3); hold on;
    p7 = plot(timex, traj.ndiff_expun_b1_CI(1, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p8 = plot(timex, traj.ndiff_expun_b1_CI(2, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p7.Color(4) = opac;
    p8.Color(4) = opac;
    
    plot(timex, zeros(1, frm), 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('Distance (a.u.)')
    n2 = title('Unexpected minus Expected')
    set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 3) %p-value
    plot(timex, traj.p_exp_b1, 'Color', diffcolor1, 'LineWidth', 2); hold on;
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    AX = subplot(3, 1, 1); %subplot 1-4: hi vs lo
    AY = subplot(3, 1, 2); %subplot 5-8: hi minus lo
    AZ = subplot(3, 1, 3); %subplot 9-12: p-values
    
    set(AX, 'YLim', [0 1.0]);
    set(AX, 'YTick', [0:0.2:1]);
    set(AX, 'XLim', [0 1500]);
    set(AX, 'XTick', [0:500:1500]);
    set(AY, 'YLim', [-0.1 0.2]);
    set(AY, 'YTick', [-0.1:0.05:0.2]);
    set(AY, 'XLim', [0 1500]);
    set(AY, 'XTick', [0:500:1500]);
    set(AZ, 'YLim', [0 0.06]);
    set(AZ, 'YTick', [0:0.02:0.06]);
    set(AZ, 'XLim', [0 1500]);
    set(AZ, 'XTick', [0:500:1500]);
end

%% 2A) resp trajectories as a function of attention
% bin 1
for nboot = 1:numboot
    fprintf('resp traj vs att_bin1: %s\n', num2str(nboot))
    sub_ind = 1:numel(allsub);
    %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
    for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
        clear nfoc_cond ndiv_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
        
        foc_cond = data.dist_foc1{sub_ind(ss)}; %all foc trials in b1 of that subj
        div_cond = data.dist_div1{sub_ind(ss)}; %all div trials in b1 of that subj
        
        %use the smaller cond (foc or div) to find index
        min_focdiv = min(size(foc_cond, 2), size(div_cond, 2));
        ind = randi(min_focdiv, [1, min_focdiv]); % this is an index matrix to sample
        
        for timepoint = 1:frm
            if size(foc_cond, 2) < size(div_cond, 2)
                sh_div = randsample(size(div_cond, 2), size(div_cond, 2))'; %shuffle trials for bootstrapping
                div_ind = sh_div(ind);
                nfoc_cond(timepoint, :) = foc_cond(timepoint, ind(1, :)); %resampled expected trials
                ndiv_cond(timepoint, :) = div_cond(timepoint, div_ind(1, :)); %resampled unexpected trials
                
            else
                sh_foc = randsample(size(foc_cond, 2), size(foc_cond, 2))'; %shuffle trials for bootstrapping
                foc_ind = sh_foc(ind);
                nfoc_cond(timepoint, :) = foc_cond(timepoint, foc_ind(1, :)); %resampled expected trials
                ndiv_cond(timepoint, :) = div_cond(timepoint, ind(1, :)); %resampled unexpected trials
            end
            nfoc_avg(timepoint, :) = nanmean(nfoc_cond(timepoint, :), 2);
            ndiv_avg(timepoint, :) = nanmean(ndiv_cond(timepoint, :), 2);
            ndiff_focdiv_avg(timepoint, :) = nfoc_avg(timepoint, :) - ndiv_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
        end
        nfoc_cond_b1{ss} = nfoc_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        ndiv_cond_b1{ss} = ndiv_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nfoc_avg_b1(:, ss) = nfoc_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        ndiv_avg_b1(:, ss) = ndiv_avg;
        ndiff_focdiv_avg_b1(:, ss) = ndiff_focdiv_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    traj.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
    traj.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
    traj.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
end

for tp = 1:frm %loop over timepoints to find average across iterations
    %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    traj.ndiff_focdiv_b1_CI(:, tp) = prctile(traj.ndiff_focdiv_savg_b1(tp, :), [2.5 97.5]);
    traj.nfoc_b1_CI(:, tp) = prctile(traj.nfoc_savg_b1(tp, :), [2.5 97.5]);
    traj.ndiv_b1_CI(:, tp) = prctile(traj.ndiv_savg_b1(tp, :), [2.5 97.5]);
    traj.p_att_b1(:, tp) = min((2*min(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot, ...
        1-(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot))), ...
        (2*min(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot, ...
        1-(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2B) plotting traj as a function of att (bin1-bin4)
for fig = 2
    figure(2); clf;
    n1 = suptitle(['Resp Traj as a f(n) of Attention (n=' num2str(size(allsub,2)) ')'])
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    %bin1
    subplot(3, 1, 1) %resampling data of exp vs unexp
    p1 = plot(timex, nanmean(traj.nfoc_savg_b1, 2), 'Color', c_foc, 'LineWidth', 3); hold on; %resampled expected trials
    p2 = plot(timex, traj.nfoc_b1_CI(1, :), 'Color', c_foc, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, traj.nfoc_b1_CI(2, :), 'Color', c_foc, 'LineWidth', 1.5); hold on;
    p2.Color(4) = opac;
    p3.Color(4) = opac;
    
    p4 = plot(timex, nanmean(traj.ndiv_savg_b1, 2), 'Color', c_div, 'LineWidth', 3); %resampled unexpected trials
    p5 = plot(timex, traj.ndiv_b1_CI(1, :), 'Color', c_div, 'LineWidth', 1.5); hold on;
    p6 = plot(timex, traj.ndiv_b1_CI(2, :), 'Color', c_div, 'LineWidth', 1.5);
    p5.Color(4) = opac;
    p6.Color(4) = opac;
    
    ylabel('Distance (a.u.)')
    legend([p1 p4], {'Focused', 'Divided'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %unexpected minus expected
    plot(timex, nanmean(traj.ndiff_focdiv_savg_b1, 2), 'Color', diffcolor1, 'LineWidth', 3); hold on;
    p7 = plot(timex, traj.ndiff_focdiv_b1_CI(1, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p8 = plot(timex, traj.ndiff_focdiv_b1_CI(2, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p7.Color(4) = opac;
    p8.Color(4) = opac;
    
    plot(timex, zeros(1, frm), 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('Distance (a.u.)')
    n2 = title('Focused minus Divided')
    set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 3) %p-value
    plot(timex, traj.p_att_b1, 'Color', diffcolor1, 'LineWidth', 2); hold on;
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    AX = subplot(3, 1, 1); %subplot 1-4: hi vs lo
    AY = subplot(3, 1, 2); %subplot 5-8: hi minus lo
    AZ = subplot(3, 1, 3); %subplot 9-12: p-values
    
    set(AX, 'YLim', [0 1.0]);
    set(AX, 'YTick', [0:0.2:1]);
    set(AX, 'XLim', [0 1500]);
    set(AX, 'XTick', [0:500:1500]);
    set(AY, 'YLim', [-0.1 0.2]);
    set(AY, 'YTick', [-0.1:0.05:0.2]);
    set(AY, 'XLim', [0 1500]);
    set(AY, 'XTick', [0:500:1500]);
    set(AZ, 'YLim', [0 0.06]);
    set(AZ, 'YTick', [0:0.02:0.06]);
    set(AZ, 'XLim', [0 1500]);
    set(AZ, 'XTick', [0:500:1500]);
end

%% 3A) resp trajectories as a function of coherence
% bin 1
for nboot = 1:numboot
    fprintf('resp traj vs coh_bin1: %s\n', num2str(nboot))
    sub_ind = 1:numel(allsub);
    %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
    for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
        clear nhi_cond nlo_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
        
        hi_cond = data.dist_hi1{sub_ind(ss)}; %all foc trials in b1 of that subj
        lo_cond = data.dist_lo1{sub_ind(ss)}; %all div trials in b1 of that subj
        
        %use the smaller cond (foc or div) to find index
        min_hilo = min(size(hi_cond, 2), size(lo_cond, 2));
        ind = randi(min_hilo, [1, min_hilo]); % this is an index matrix to sample
        
        for timepoint = 1:frm
            if size(hi_cond, 2) < size(lo_cond, 2)
                sh_lo = randsample(size(lo_cond, 2), size(lo_cond, 2))'; %shuffle trials for bootstrapping
                lo_ind = sh_lo(ind);
                nhi_cond(timepoint, :) = hi_cond(timepoint, ind(1, :)); %resampled expected trials
                nlo_cond(timepoint, :) = lo_cond(timepoint, lo_ind(1, :)); %resampled unexpected trials
                
            else
                sh_hi = randsample(size(hi_cond, 2), size(hi_cond, 2))'; %shuffle trials for bootstrapping
                hi_ind = sh_hi(ind);
                nhi_cond(timepoint, :) = hi_cond(timepoint, hi_ind(1, :)); %resampled expected trials
                nlo_cond(timepoint, :) = lo_cond(timepoint, ind(1, :)); %resampled unexpected trials
            end
            nhi_avg(timepoint, :) = nanmean(nhi_cond(timepoint, :), 2);
            nlo_avg(timepoint, :) = nanmean(nlo_cond(timepoint, :), 2);
            ndiff_hilo_avg(timepoint, :) = nhi_avg(timepoint, :) - nlo_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
        end
        nhi_cond_b1{ss} = nhi_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        nlo_cond_b1{ss} = nlo_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nhi_avg_b1(:, ss) = nhi_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        nlo_avg_b1(:, ss) = nlo_avg;
        ndiff_hilo_avg_b1(:, ss) = ndiff_hilo_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    traj.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
    traj.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
    traj.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
end

for tp = 1:frm %loop over timepoints to find average across iterations
    %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    traj.ndiff_hilo_b1_CI(:, tp) = prctile(traj.ndiff_hilo_savg_b1(tp, :), [2.5 97.5]);
    traj.nhi_b1_CI(:, tp) = prctile(traj.nhi_savg_b1(tp, :), [2.5 97.5]);
    traj.nlo_b1_CI(:, tp) = prctile(traj.nlo_savg_b1(tp, :), [2.5 97.5]);
    traj.p_coh_b1(:, tp) = min((2*min(numel(find(traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot, ...
        1-(numel(find(traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot))), ...
        (2*min(numel(find(traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot, ...
        1-(numel(find(traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot))));
end

%% 3B) plotting traj as a function of coh (bin1-bin4)
for fig = 3
    figure(3); clf;
    n1 = suptitle(['Resp Traj as a f(n) of Coherence (n=' num2str(size(allsub,2)) ')'])
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    %bin1
    subplot(3, 1, 1) %resampling data of hi vs lo
    p1 = plot(timex, nanmean(traj.nhi_savg_b1, 2), 'Color', c_hi, 'LineWidth', 3); hold on; %resampled expected trials
    p2 = plot(timex, traj.nhi_b1_CI(1, :), 'Color', c_hi, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, traj.nhi_b1_CI(2, :), 'Color', c_hi, 'LineWidth', 1.5); hold on;
    p2.Color(4) = opac;
    p3.Color(4) = opac;
    
    p4 = plot(timex, nanmean(traj.nlo_savg_b1, 2), 'Color', c_lo, 'LineWidth', 3); %resampled unexpected trials
    p5 = plot(timex, traj.nlo_b1_CI(1, :), 'Color', c_lo, 'LineWidth', 1.5); hold on;
    p6 = plot(timex, traj.nlo_b1_CI(2, :), 'Color', c_lo, 'LineWidth', 1.5);
    p5.Color(4) = opac;
    p6.Color(4) = opac;
    
    ylabel('Distance (a.u.)')
    legend([p1 p4], {'Hi Coh', 'Lo Coh'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %hi minus lo
    plot(timex, nanmean(traj.ndiff_hilo_savg_b1, 2), 'Color', diffcolor1, 'LineWidth', 3); hold on;
    p7 = plot(timex, traj.ndiff_hilo_b1_CI(1, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p8 = plot(timex, traj.ndiff_hilo_b1_CI(2, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p7.Color(4) = opac;
    p8.Color(4) = opac;
    
    plot(timex, zeros(1, frm), 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('Distance (a.u.)')
    n2 = title('Hi Coh minus Lo Coh')
    set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 3) %p-value
    plot(timex, traj.p_coh_b1, 'Color', diffcolor1, 'LineWidth', 2); hold on;
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    
    AX = subplot(3, 1, 1); %subplot 1-4: hi vs lo
    AY = subplot(3, 1, 2); %subplot 5-8: hi minus lo
    AZ = subplot(3, 1, 3); %subplot 9-12: p-values
    
    set(AX, 'YLim', [0 1.0]);
    set(AX, 'YTick', [0:0.2:1]);
    set(AX, 'XLim', [0 1500]);
    set(AX, 'XTick', [0:500:1500]);
    set(AY, 'YLim', [-0.1 0.2]);
    set(AY, 'YTick', [-0.1:0.05:0.2]);
    set(AY, 'XLim', [0 1500]);
    set(AY, 'XTick', [0:500:1500]);
    set(AZ, 'YLim', [0 0.06]);
    set(AZ, 'YTick', [0:0.02:0.06]);
    set(AZ, 'XLim', [0 1500]);
    set(AZ, 'XTick', [0:500:1500]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%% Resp Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4A) resp error as a function of expectation
% bin 1
for nboot = 1:numboot
    fprintf('resp error vs exp_bin1: %s\n', num2str(nboot))
    sub_ind = 1:numel(allsub);
    %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
    for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
        clear nexp_cond nun_cond ndiff_expun_avg %doing this b/c #of trials in each bin vary across subj
        
        exp_cond = data.resp_exp1{sub_ind(ss)}; %all exp trials in b1 of that subj
        un_cond = data.resp_un1{sub_ind(ss)}; %all un trials in b1 of that subj
        
        %%% 2 ways of doing this i think: 1) fixed # across subj determined by the
        %%% min # of un trials across subj; 2) the number of un of each subj
        %%% trying 2) first
        ind = randi(size(un_cond, 2), [1, size(un_cond, 2)]); % this is an index matrix to sample
        % out XX trials (with replacement) each from expected and
        % unexpected condition
        for timepoint = 1:frm
            sh_exp = randsample(size(exp_cond, 2), size(exp_cond, 2))'; %shuffle trials for bootstrapping
            exp_ind = sh_exp(ind);
            nexp_cond(timepoint, :) = exp_cond(timepoint, exp_ind(1, :)); %resampled expected trials
            nexp_avg(timepoint, :) = nanmean(nexp_cond(timepoint, :), 2);
            nun_cond(timepoint, :) = un_cond(timepoint, ind(1, :)); %resampled unexpected trials
            nun_avg(timepoint, :) = nanmean(nun_cond(timepoint, :), 2);
            ndiff_expun_avg(timepoint, :) = nun_avg(timepoint, :) - nexp_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
        end
        nexp_cond_b1{ss} = nexp_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        nun_cond_b1{ss} = nun_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nexp_avg_b1(:, ss) = nexp_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        nun_avg_b1(:, ss) = nun_avg;
        ndiff_expun_avg_b1(:, ss) = ndiff_expun_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    resp.nexp_savg_b1(:, nboot) = nanmean(nexp_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
    resp.nun_savg_b1(:, nboot) = nanmean(nun_avg_b1, 2);
    resp.ndiff_expun_savg_b1(:, nboot) = nanmean(ndiff_expun_avg_b1, 2);
end

for tp = 1:frm %loop over timepoints to find average across iterations
    %display([nanmean(resp.ndiff_savg_b1(tp, :)), prctile(resp.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    resp.ndiff_expun_b1_CI(:, tp) = prctile(resp.ndiff_expun_savg_b1(tp, :), [2.5 97.5]);
    resp.nexp_b1_CI(:, tp) = prctile(resp.nexp_savg_b1(tp, :), [2.5 97.5]);
    resp.nun_b1_CI(:, tp) = prctile(resp.nun_savg_b1(tp, :), [2.5 97.5]);
    resp.p_exp_b1(:, tp) = min((2*min(numel(find(resp.ndiff_expun_savg_b1(tp, :) > 0))./numboot, ...
        1-(numel(find(resp.ndiff_expun_savg_b1(tp, :) < 0))./numboot))), ...
        (2*min(numel(find(resp.ndiff_expun_savg_b1(tp, :) < 0))./numboot, ...
        1-(numel(find(resp.ndiff_expun_savg_b1(tp, :) > 0))./numboot))));
    
    
end

%% 4B) plotting resp error as a function of exp (bin1-bin4)
for fig = 4
    figure(4); clf;
    n1 = suptitle(['Resp Error as a f(n) of Expectation (n=' num2str(size(allsub,2)) ')'])
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    %bin1
    subplot(3, 1, 1) %resampling data of exp vs unexp
    p1 = plot(timex, nanmean(resp.nexp_savg_b1, 2), 'Color', c_exp, 'LineWidth', 3); hold on; %resampled expected trials
    p2 = plot(timex, resp.nexp_b1_CI(1, :), 'Color', c_exp, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, resp.nexp_b1_CI(2, :), 'Color', c_exp, 'LineWidth', 1.5); hold on;
    p2.Color(4) = opac;
    p3.Color(4) = opac;
    
    p4 = plot(timex, nanmean(resp.nun_savg_b1, 2), 'Color', c_unexp, 'LineWidth', 3); %resampled unexpected trials
    p5 = plot(timex, resp.nun_b1_CI(1, :), 'Color', c_unexp, 'LineWidth', 1.5); hold on;
    p6 = plot(timex, resp.nun_b1_CI(2, :), 'Color', c_unexp, 'LineWidth', 1.5);
    p5.Color(4) = opac;
    p6.Color(4) = opac;
    
    ylabel('Degree Error (deg)')
    legend([p1 p4], {'Expected', 'Unexpected'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %unexpected minus expected
    plot(timex, nanmean(resp.ndiff_expun_savg_b1, 2), 'Color', diffcolor1, 'LineWidth', 3); hold on;
    p7 = plot(timex, resp.ndiff_expun_b1_CI(1, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p8 = plot(timex, resp.ndiff_expun_b1_CI(2, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p7.Color(4) = opac;
    p8.Color(4) = opac;
    
    plot(timex, zeros(1, frm), 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('Degree Error (deg)')
    n2 = title('Unexpected minus Expected')
    set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 3) %p-value
    plot(timex, resp.p_exp_b1, 'Color', diffcolor1, 'LineWidth', 2); hold on;
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    AX = subplot(3, 1, 1); %subplot 1-4: hi vs lo
    AY = subplot(3, 1, 2); %subplot 5-8: hi minus lo
    AZ = subplot(3, 1, 3); %subplot 9-12: p-values
    
    set(AX, 'YLim', [0 100]);
    set(AX, 'YTick', [0:20:100]);
    set(AX, 'XLim', [0 1500]);
    set(AX, 'XTick', [0:500:1500]);
    set(AY, 'YLim', [-10 30]);
    set(AY, 'YTick', [-10:10:30]);
    set(AY, 'XLim', [0 1500]);
    set(AY, 'XTick', [0:500:1500]);
    set(AZ, 'YLim', [0 0.06]);
    set(AZ, 'YTick', [0:0.02:0.06]);
    set(AZ, 'XLim', [0 1500]);
    set(AZ, 'XTick', [0:500:1500]);
end

%% 5A) resp error as a function of attention
% bin 1
for nboot = 1:numboot
    fprintf('resp error vs att_bin1: %s\n', num2str(nboot))
    sub_ind = 1:numel(allsub);
    %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
    for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
        clear nfoc_cond ndiv_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
        
        foc_cond = data.resp_foc1{sub_ind(ss)}; %all foc trials in b1 of that subj
        div_cond = data.resp_div1{sub_ind(ss)}; %all div trials in b1 of that subj
        
        %use the smaller cond (foc or div) to find index
        min_focdiv = min(size(foc_cond, 2), size(div_cond, 2));
        ind = randi(min_focdiv, [1, min_focdiv]); % this is an index matrix to sample
        
        for timepoint = 1:frm
            if size(foc_cond, 2) < size(div_cond, 2)
                sh_div = randsample(size(div_cond, 2), size(div_cond, 2))'; %shuffle trials for bootstrapping
                div_ind = sh_div(ind);
                nfoc_cond(timepoint, :) = foc_cond(timepoint, ind(1, :)); %resampled expected trials
                ndiv_cond(timepoint, :) = div_cond(timepoint, div_ind(1, :)); %resampled unexpected trials
                
            else
                sh_foc = randsample(size(foc_cond, 2), size(foc_cond, 2))'; %shuffle trials for bootstrapping
                foc_ind = sh_foc(ind);
                nfoc_cond(timepoint, :) = foc_cond(timepoint, foc_ind(1, :)); %resampled expected trials
                ndiv_cond(timepoint, :) = div_cond(timepoint, ind(1, :)); %resampled unexpected trials
            end
            nfoc_avg(timepoint, :) = nanmean(nfoc_cond(timepoint, :), 2);
            ndiv_avg(timepoint, :) = nanmean(ndiv_cond(timepoint, :), 2);
            ndiff_focdiv_avg(timepoint, :) = ndiv_avg(timepoint, :) - nfoc_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
        end
        nfoc_cond_b1{ss} = nfoc_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        ndiv_cond_b1{ss} = ndiv_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nfoc_avg_b1(:, ss) = nfoc_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        ndiv_avg_b1(:, ss) = ndiv_avg;
        ndiff_focdiv_avg_b1(:, ss) = ndiff_focdiv_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    resp.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
    resp.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
    resp.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
end

for tp = 1:frm %loop over timepoints to find average across iterations
    %display([nanmean(resp.ndiff_savg_b1(tp, :)), prctile(resp.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    resp.ndiff_focdiv_b1_CI(:, tp) = prctile(resp.ndiff_focdiv_savg_b1(tp, :), [2.5 97.5]);
    resp.nfoc_b1_CI(:, tp) = prctile(resp.nfoc_savg_b1(tp, :), [2.5 97.5]);
    resp.ndiv_b1_CI(:, tp) = prctile(resp.ndiv_savg_b1(tp, :), [2.5 97.5]);
    resp.p_att_b1(:, tp) = min((2*min(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot, ...
        1-(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot))), ...
        (2*min(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot, ...
        1-(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot))));
    
end

%% 5B) plotting resp error as a function of att (bin1-bin4)
for fig = 5
    figure(5); clf;
    n1 = suptitle(['Resp Error as a f(n) of Attention (n=' num2str(size(allsub,2)) ')'])
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    %bin1
    subplot(3, 1, 1) %resampling data of exp vs unexp
    p1 = plot(timex, nanmean(resp.nfoc_savg_b1, 2), 'Color', c_foc, 'LineWidth', 3); hold on; %resampled expected trials
    p2 = plot(timex, resp.nfoc_b1_CI(1, :), 'Color', c_foc, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, resp.nfoc_b1_CI(2, :), 'Color', c_foc, 'LineWidth', 1.5); hold on;
    p2.Color(4) = opac;
    p3.Color(4) = opac;
    
    p4 = plot(timex, nanmean(resp.ndiv_savg_b1, 2), 'Color', c_div, 'LineWidth', 3); %resampled unexpected trials
    p5 = plot(timex, resp.ndiv_b1_CI(1, :), 'Color', c_div, 'LineWidth', 1.5); hold on;
    p6 = plot(timex, resp.ndiv_b1_CI(2, :), 'Color', c_div, 'LineWidth', 1.5);
    p5.Color(4) = opac;
    p6.Color(4) = opac;
    
    ylabel('Degree Error (deg)')
    legend([p1 p4], {'Focused', 'Divided'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %unexpected minus expected
    plot(timex, nanmean(resp.ndiff_focdiv_savg_b1, 2), 'Color', diffcolor1, 'LineWidth', 3); hold on;
    p7 = plot(timex, resp.ndiff_focdiv_b1_CI(1, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p8 = plot(timex, resp.ndiff_focdiv_b1_CI(2, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p7.Color(4) = opac;
    p8.Color(4) = opac;
    
    plot(timex, zeros(1, frm), 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('Degree Error (deg)')
    n2 = title('Focused minus Divided')
    set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 3) %p-value
    plot(timex, resp.p_att_b1, 'Color', diffcolor1, 'LineWidth', 2); hold on;
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    AX = subplot(3, 1, 1); %subplot 1-4: hi vs lo
    AY = subplot(3, 1, 2); %subplot 5-8: hi minus lo
    AZ = subplot(3, 1, 3); %subplot 9-12: p-values
    
    set(AX, 'YLim', [0 100]);
    set(AX, 'YTick', [0:20:100]);
    set(AX, 'XLim', [0 1500]);
    set(AX, 'XTick', [0:500:1500]);
    set(AY, 'YLim', [-10 30]);
    set(AY, 'YTick', [-10:10:30]);
    set(AY, 'XLim', [0 1500]);
    set(AY, 'XTick', [0:500:1500]);
    set(AZ, 'YLim', [0 0.06]);
    set(AZ, 'YTick', [0:0.02:0.06]);
    set(AZ, 'XLim', [0 1500]);
    set(AZ, 'XTick', [0:500:1500]);
end

%% 6A) resp trajectories as a function of coherence
% bin 1
for nboot = 1:numboot
    fprintf('resp error vs coh_bin1: %s\n', num2str(nboot))
    sub_ind = 1:numel(allsub);
    %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
    for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
        clear nhi_cond nlo_cond ndiff_hilo_avg %doing this b/c #of trials in each bin vary across subj
        
        hi_cond = data.resp_hi1{sub_ind(ss)}; %all foc trials in b1 of that subj
        lo_cond = data.resp_lo1{sub_ind(ss)}; %all div trials in b1 of that subj
        
        %use the smaller cond (foc or div) to find index
        min_hilo = min(size(hi_cond, 2), size(lo_cond, 2));
        ind = randi(min_hilo, [1, min_hilo]); % this is an index matrix to sample
        
        for timepoint = 1:frm
            if size(hi_cond, 2) < size(lo_cond, 2)
                sh_lo = randsample(size(lo_cond, 2), size(lo_cond, 2))'; %shuffle trials for bootstrapping
                lo_ind = sh_lo(ind);
                nhi_cond(timepoint, :) = hi_cond(timepoint, ind(1, :)); %resampled expected trials
                nlo_cond(timepoint, :) = lo_cond(timepoint, lo_ind(1, :)); %resampled unexpected trials
                
            else
                sh_hi = randsample(size(hi_cond, 2), size(hi_cond, 2))'; %shuffle trials for bootstrapping
                hi_ind = sh_hi(ind);
                nhi_cond(timepoint, :) = hi_cond(timepoint, hi_ind(1, :)); %resampled expected trials
                nlo_cond(timepoint, :) = lo_cond(timepoint, ind(1, :)); %resampled unexpected trials
            end
            nhi_avg(timepoint, :) = nanmean(nhi_cond(timepoint, :), 2);
            nlo_avg(timepoint, :) = nanmean(nlo_cond(timepoint, :), 2);
            %ndiff_hilo_avg(timepoint, :) = nhi_avg(timepoint, :) - nlo_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
            ndiff_hilo_avg(timepoint, :) = nlo_avg(timepoint, :) - nhi_avg(timepoint, :);
        end
        nhi_cond_b1{ss} = nhi_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        nlo_cond_b1{ss} = nlo_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nhi_avg_b1(:, ss) = nhi_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        nlo_avg_b1(:, ss) = nlo_avg;
        ndiff_hilo_avg_b1(:, ss) = ndiff_hilo_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    resp.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
    resp.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
    resp.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
end

for tp = 1:frm %loop over timepoints to find average across iterations
    %display([nanmean(resp.ndiff_savg_b1(tp, :)), prctile(resp.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    resp.ndiff_hilo_b1_CI(:, tp) = prctile(resp.ndiff_hilo_savg_b1(tp, :), [2.5 97.5]);
    resp.nhi_b1_CI(:, tp) = prctile(resp.nhi_savg_b1(tp, :), [2.5 97.5]);
    resp.nlo_b1_CI(:, tp) = prctile(resp.nlo_savg_b1(tp, :), [2.5 97.5]);
    resp.p_coh_b1(:, tp) = min((2*min(numel(find(resp.ndiff_hilo_savg_b1(tp, :) > 0))./numboot, ...
        1-(numel(find(resp.ndiff_hilo_savg_b1(tp, :) < 0))./numboot))), ...
        (2*min(numel(find(resp.ndiff_hilo_savg_b1(tp, :) < 0))./numboot, ...
        1-(numel(find(resp.ndiff_hilo_savg_b1(tp, :) > 0))./numboot))));
end

%% 6B) plotting traj as a function of coh (bin1-bin4)
for fig = 6
    figure(6); clf;
    n1 = suptitle(['Resp Error as a f(n) of Coherence (n=' num2str(size(allsub,2)) ')'])
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    %bin1
    subplot(3, 1, 1) %resampling data of hi vs lo
    p1 = plot(timex, nanmean(resp.nhi_savg_b1, 2), 'Color', c_hi, 'LineWidth', 3); hold on; %resampled expected trials
    p2 = plot(timex, resp.nhi_b1_CI(1, :), 'Color', c_hi, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, resp.nhi_b1_CI(2, :), 'Color', c_hi, 'LineWidth', 1.5); hold on;
    p2.Color(4) = opac;
    p3.Color(4) = opac;
    
    p4 = plot(timex, nanmean(resp.nlo_savg_b1, 2), 'Color', c_lo, 'LineWidth', 3); %resampled unexpected trials
    p5 = plot(timex, resp.nlo_b1_CI(1, :), 'Color', c_lo, 'LineWidth', 1.5); hold on;
    p6 = plot(timex, resp.nlo_b1_CI(2, :), 'Color', c_lo, 'LineWidth', 1.5);
    p5.Color(4) = opac;
    p6.Color(4) = opac;
    
    ylabel('Degree Error (deg)')
    legend([p1 p4], {'Hi Coh', 'Lo Coh'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %hi minus lo
    plot(timex, nanmean(resp.ndiff_hilo_savg_b1, 2), 'Color', diffcolor1, 'LineWidth', 3); hold on;
    p7 = plot(timex, resp.ndiff_hilo_b1_CI(1, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p8 = plot(timex, resp.ndiff_hilo_b1_CI(2, :), 'Color', diffcolor1, 'LineWidth', 1.5); hold on;
    p7.Color(4) = opac;
    p8.Color(4) = opac;
    
    plot(timex, zeros(1, frm), 'k--', 'LineWidth', 2);
    
    xlabel('Time (ms)')
    ylabel('Degree Error (deg)')
    n2 = title('Hi Coh minus Lo Coh')
    set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 3) %p-value
    plot(timex, resp.p_coh_b1, 'Color', diffcolor1, 'LineWidth', 2); hold on;
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    AX = subplot(3, 1, 1); %subplot 1-4: hi vs lo
    AY = subplot(3, 1, 2); %subplot 5-8: hi minus lo
    AZ = subplot(3, 1, 3); %subplot 9-12: p-values
    
    set(AX, 'YLim', [0 100]);
    set(AX, 'YTick', [0:20:100]);
    set(AX, 'XLim', [0 1500]);
    set(AX, 'XTick', [0:500:1500]);
    set(AY, 'YLim', [-10 30]);
    set(AY, 'YTick', [-10:10:30]);
    set(AY, 'XLim', [0 1500]);
    set(AY, 'XTick', [0:500:1500]);
    set(AZ, 'YLim', [0 0.06]);
    set(AZ, 'YTick', [0:0.02:0.06]);
    set(AZ, 'XLim', [0 1500]);
    set(AZ, 'XTick', [0:500:1500]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%save('bootstrap_controls.mat', 'traj', 'resp')
save('bootstrap_patients.mat', 'traj', 'resp')

toc