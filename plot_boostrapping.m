% This script plots the results

close all; clear all;
load('controls_data.mat')
load('bootstrap_controls.mat')

ctraj = traj;
cresp = resp;
cdata = data;
callsub = allsub;
clear traj resp data allsub

load('patients_data.mat')
load('bootstrap_patients.mat')

ptraj = traj;
presp = resp;
pdata = data;
pallsub = allsub;

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

c1 = [0.2 0.2 0.2]; %dark grey: coh
c2 = [0.2 0.4 1]; %blue: att
c3 = [0.8 0 0]; %red: exp

%response trajectories
for fig = 1
    figure(1); clf;
    n1 = suptitle(['Response Trajectories (10 controls (left); 4 patients (right))'])
    %n1 = suptitle(['Controls: Response Trajectories (n=' num2str(size(allsub,2)) ')'])
    %n1 = suptitle(['Patients: Response Trajectories (n=' num2str(size(allsub,2)) ')'])
    
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    subplot(5, 2, 1) %control_coh_bin1
    p1 = plot(timex, nanmean(ctraj.nhi_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(ctraj.nlo_savg_b1, 2), 'Color', c1, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    %ylabel('Distance (a.u.)')
    legend([p1 p4], {'Hi coh', 'Lo coh'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 2) %patients_coh_bin1
    p1 = plot(timex, nanmean(ptraj.nhi_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(ptraj.nlo_savg_b1, 2), 'Color', c1, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    %ylabel('Distance (a.u.)')
    legend([p1 p4], {'Hi coh', 'Lo coh'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 3) %control_att_bin1
    p1 = plot(timex, nanmean(ctraj.nfoc_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(ctraj.ndiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    ylabel('Distance (a.u.)')
    legend([p1 p4], {'Focused', 'Divided'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    %
    subplot(5, 2, 4) %pt_att_bin1
    p1 = plot(timex, nanmean(ptraj.nfoc_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(ptraj.ndiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    ylabel('Distance (a.u.)')
    legend([p1 p4], {'Focused', 'Divided'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 5) %control_exp_bin1
    p1 = plot(timex, nanmean(ctraj.nexp_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(ctraj.nun_savg_b1, 2), 'Color', c3, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    legend([p1 p4], {'Expected', 'Unexpected'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 6) %patients_exp_bin1
    p1 = plot(timex, nanmean(ptraj.nexp_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(ptraj.nun_savg_b1, 2), 'Color', c3, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    legend([p1 p4], {'Expected', 'Unexpected'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 7) %control_diff_bin1
    %hi coh-lo coh
    m1 = plot(timex, nanmean(ctraj.ndiff_hilo_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on;
    m2 = plot(timex, ctraj.ndiff_hilo_b1_CI(1, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m3 = plot(timex, ctraj.ndiff_hilo_b1_CI(2, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m2.Color(4) = opac;
    m3.Color(4) = opac;
    ylabel('Distance (a.u.)')
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %foc-div
    m4 = plot(timex, nanmean(ctraj.ndiff_focdiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on;
    m5 = plot(timex, ctraj.ndiff_focdiv_b1_CI(1, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m6 = plot(timex, ctraj.ndiff_focdiv_b1_CI(2, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m5.Color(4) = opac;
    m6.Color(4) = opac;
    
    %exp-unexp
    m7 = plot(timex, (-1)*nanmean(ctraj.ndiff_expun_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on;
    m8 = plot(timex, (-1)*ctraj.ndiff_expun_b1_CI(1, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m9 = plot(timex, (-1)*ctraj.ndiff_expun_b1_CI(2, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m8.Color(4) = opac;
    m9.Color(4) = opac;
    %
    subplot(5, 2, 8) %patients_diff_bin1
    %hi coh-lo coh
    m1 = plot(timex, nanmean(ptraj.ndiff_hilo_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on;
    m2 = plot(timex, ptraj.ndiff_hilo_b1_CI(1, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m3 = plot(timex, ptraj.ndiff_hilo_b1_CI(2, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m2.Color(4) = opac;
    m3.Color(4) = opac;
    ylabel('Distance (a.u.)')
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %foc-div
    m4 = plot(timex, nanmean(ptraj.ndiff_focdiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on;
    m5 = plot(timex, ptraj.ndiff_focdiv_b1_CI(1, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m6 = plot(timex, ptraj.ndiff_focdiv_b1_CI(2, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m5.Color(4) = opac;
    m6.Color(4) = opac;
    
    %exp-unexp
    m7 = plot(timex, (-1)*nanmean(ptraj.ndiff_expun_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on;
    m8 = plot(timex, (-1)*ptraj.ndiff_expun_b1_CI(1, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m9 = plot(timex, (-1)*ptraj.ndiff_expun_b1_CI(2, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m8.Color(4) = opac;
    m9.Color(4) = opac;
    
    subplot(5, 2, 9) %control_p_bin1
    plot(timex, ctraj.p_coh_b1, 'Color', c1, 'LineWidth', 2); hold on;
    plot(timex, ctraj.p_att_b1, 'Color', c2, 'LineWidth', 2); hold on;
    plot(timex, ctraj.p_exp_b1, 'Color', c3, 'LineWidth', 2); hold on;
    
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    legend([m1 m4 m7], {'Hi-Lo', 'Foc-Div', 'Exp-Unexp'});
    ylim([0 0.06])
    xlim([0 1500])
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 10) %patients_p_bin1
    plot(timex, ptraj.p_coh_b1, 'Color', c1, 'LineWidth', 2); hold on;
    plot(timex, ptraj.p_att_b1, 'Color', c2, 'LineWidth', 2); hold on;
    plot(timex, ptraj.p_exp_b1, 'Color', c3, 'LineWidth', 2); hold on;
    
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    legend([m1 m4 m7], {'Hi-Lo', 'Foc-Div', 'Exp-Unexp'});
    ylim([0 0.06])
    xlim([0 1500])
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    for n = 1:2
        AM(n) = subplot(5, 2, n); %hi vs low
        AN(n) = subplot(5, 2, n+2); %foc vs div
        AO(n) = subplot(5, 2, n+4); %exp vs unexp
        AP(n) = subplot(5, 2, n+6); %diff plots
        AQ(n) = subplot(5, 2, n+8) %p-values plots
        
    end
    set([AM AN AO], 'YLim', [0 0.9]);
    set([AM AN AO], 'YTick', [0:0.2:0.9]);
    set([AM AN AO], 'XLim', [0 1500]);
    set([AM AN AO], 'XTick', [0:500:1500]);
    set(AP, 'YLim', [-0.2 0.2]);
    set(AP, 'YTick', [-0.2:0.05:0.2]);
    set(AP, 'XLim', [0 1500]);
    set(AP, 'XTick', [0:500:1500]);
    set(AQ, 'YLim', [0 0.06]);
    set(AQ, 'YTick', [0:0.02:0.06]);
    set(AQ, 'XLim', [0 1500]);
    set(AQ, 'XTick', [0:500:1500]);
end

%response error
for fig = 2
    figure(2); clf;
    n1 = suptitle(['Response Errors (10 controls (left); 4 patients (right))'])
    %n1 = suptitle(['Controls: Response Error (n=' num2str(size(allsub,2)) ')'])
    %n1 = suptitle(['Patients: Response Error (n=' num2str(size(allsub,2)) ')'])
    
    set(n1, 'Fontsize', 13, 'fontWeight', 'bold')
    
    subplot(5, 2, 1) %control:coh_bin1
    p1 = plot(timex, nanmean(cresp.nhi_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(cresp.nlo_savg_b1, 2), 'Color', c1, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    %ylabel('Degree Error (deg)')
    legend([p1 p4], {'Hi Coh', 'Lo Coh'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    %
    subplot(5, 2, 2) %patient:coh_bin1
    p1 = plot(timex, nanmean(presp.nhi_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(presp.nlo_savg_b1, 2), 'Color', c1, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    %ylabel('Degree Error (deg)')
    legend([p1 p4], {'Hi Coh', 'Lo Coh'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    %
    subplot(5, 2, 3) %control_att_bin1
    p1 = plot(timex, nanmean(cresp.nfoc_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(cresp.ndiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    ylabel('Degree Error (deg)')
    legend([p1 p4], {'Focused', 'Divided'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 4) %patients_att_bin1
    p1 = plot(timex, nanmean(presp.nfoc_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(presp.ndiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    ylabel('Degree Error (deg)')
    legend([p1 p4], {'Focused', 'Divided'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    %
    subplot(5, 2, 5) %control_exp_bin1
    p1 = plot(timex, nanmean(cresp.nexp_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(cresp.nun_savg_b1, 2), 'Color', c3, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    %ylabel('Degree Error (deg)')
    legend([p1 p4], {'Expected', 'Unexpected'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 6) %patients_exp_bin1
    p1 = plot(timex, nanmean(presp.nexp_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on; %resampled expected trials
    p4 = plot(timex, nanmean(presp.nun_savg_b1, 2), 'Color', c3, 'LineWidth', 2, 'LineStyle', ':'); %resampled unexpected trials
    %ylabel('Degree Error (deg)')
    legend([p1 p4], {'Expected', 'Unexpected'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(5, 2, 7) %diff_bin1
    %hi coh-lo coh
    m1 = plot(timex, (-1)*nanmean(cresp.ndiff_hilo_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on;
    m2 = plot(timex, (-1)*cresp.ndiff_hilo_b1_CI(1, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m3 = plot(timex, (-1)*cresp.ndiff_hilo_b1_CI(2, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m2.Color(4) = opac;
    m3.Color(4) = opac;
    ylabel('Degree Error (deg)')
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %foc-div
    m4 = plot(timex, (-1)*nanmean(cresp.ndiff_focdiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on;
    m5 = plot(timex, (-1)*cresp.ndiff_focdiv_b1_CI(1, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m6 = plot(timex, (-1)*cresp.ndiff_focdiv_b1_CI(2, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m5.Color(4) = opac;
    m6.Color(4) = opac;
    
    %exp-unexp
    m7 = plot(timex, (-1)*nanmean(cresp.ndiff_expun_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on;
    m8 = plot(timex, (-1)*cresp.ndiff_expun_b1_CI(1, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m9 = plot(timex, (-1)*cresp.ndiff_expun_b1_CI(2, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m8.Color(4) = opac;
    m9.Color(4) = opac;
    
    subplot(5, 2, 8) %diff_bin1
    %hi coh-lo coh
    m1 = plot(timex, (-1)*nanmean(presp.ndiff_hilo_savg_b1, 2), 'Color', c1, 'LineWidth', 2); hold on;
    m2 = plot(timex, (-1)*presp.ndiff_hilo_b1_CI(1, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m3 = plot(timex, (-1)*presp.ndiff_hilo_b1_CI(2, :), 'Color', c1, 'LineWidth', 1.5); hold on;
    m2.Color(4) = opac;
    m3.Color(4) = opac;
    ylabel('Degree Error (deg)')
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %foc-div
    m4 = plot(timex, (-1)*nanmean(presp.ndiff_focdiv_savg_b1, 2), 'Color', c2, 'LineWidth', 2); hold on;
    m5 = plot(timex, (-1)*presp.ndiff_focdiv_b1_CI(1, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m6 = plot(timex, (-1)*presp.ndiff_focdiv_b1_CI(2, :), 'Color', c2, 'LineWidth', 1.5); hold on;
    m5.Color(4) = opac;
    m6.Color(4) = opac;
    
    %exp-unexp
    m7 = plot(timex, (-1)*nanmean(presp.ndiff_expun_savg_b1, 2), 'Color', c3, 'LineWidth', 2); hold on;
    m8 = plot(timex, (-1)*presp.ndiff_expun_b1_CI(1, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m9 = plot(timex, (-1)*presp.ndiff_expun_b1_CI(2, :), 'Color', c3, 'LineWidth', 1.5); hold on;
    m8.Color(4) = opac;
    m9.Color(4) = opac;
    
    subplot(5, 2, 9) %control_p_bin1
    plot(timex, cresp.p_coh_b1, 'Color', c1, 'LineWidth', 2); hold on;
    plot(timex, cresp.p_att_b1, 'Color', c2, 'LineWidth', 2); hold on;
    plot(timex, cresp.p_exp_b1, 'Color', c3, 'LineWidth', 2); hold on;
    
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    legend([m1 m4 m7], {'Hi-Lo', 'Foc-Div', 'Exp-Unexp'});
    ylim([0 0.06])
    xlim([0 1500])
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    
    subplot(5, 2, 10) %p_bin1
    plot(timex, presp.p_coh_b1, 'Color', c1, 'LineWidth', 2); hold on;
    plot(timex, presp.p_att_b1, 'Color', c2, 'LineWidth', 2); hold on;
    plot(timex, presp.p_exp_b1, 'Color', c3, 'LineWidth', 2); hold on;
    
    plot(timex, ones(1, frm).*0.05, 'k--', 'LineWidth', 2);
    legend([m1 m4 m7], {'Hi-Lo', 'Foc-Div', 'Exp-Unexp'});
    ylim([0 0.06])
    xlim([0 1500])
    xlabel('Time (ms)')
    ylabel('p-value')
    set(n1, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    for n = 1:2
        AM(n) = subplot(5, 2, n); %hi vs low
        AN(n) = subplot(5, 2, n+2); %foc vs div
        AO(n) = subplot(5, 2, n+4); %exp vs unexp
        AP(n) = subplot(5, 2, n+6); %diff plots
        AQ(n) = subplot(5, 2, n+8) %p-values plots
    end
    
    %         set([AM AN AO], 'YLim', [0 100]);
    %         set([AM AN AO], 'YTick', [0:20:100]);
    set([AM AN AO], 'YLim', [20 90]);
    set([AM AN AO], 'YTick', [20:20:90])
    set([AM AN AO], 'XLim', [0 1500]);
    set([AM AN AO], 'XTick', [0:500:1500]);
    set(AP, 'YLim', [-20 10]);
    set(AP, 'YTick', [-20:10:10]);
    set(AP, 'XLim', [0 1500]);
    set(AP, 'XTick', [0:500:1500]);
    set(AQ, 'YLim', [0 0.06]);
    set(AQ, 'YTick', [0:0.02:0.06]);
    set(AQ, 'XLim', [0 1500]);
    set(AQ, 'XTick', [0:500:1500]);
end

