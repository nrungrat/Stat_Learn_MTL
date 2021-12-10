% This script analyzes individual-subject response errors

close all; clear all;
tic

%for comparison:
load('controls_data_mat')
load('bootstrap_controls.mat')

ctraj = traj;
cresp = resp;
cdata = data;
callsub = allsub;
clear traj resp data allsub

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
c_pat = [0 0.4 0.4]; %ocean green
c_con = [0.3137 0.3137 0.3137]; %dark grey
opac = 0.6; %opacity value for plotting CI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss = 1:size(allsub, 2)
    ss
    clear traj resp
    %%%%%%%%%%%%%%%%%%%resp_expectation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nboot = 1:numboot
        clear nexp_cond nun_cond ndiff_expun_avg
        
        exp_cond = data.resp_exp1{ss}; %all exp trials in b1 of that subj
        un_cond = data.resp_un1{ss}; %all un trials in b1 of that subj
        ind = randi(size(un_cond, 2), [1, size(un_cond, 2)]); % this is an index matrix to sample
        for timepoint = 1:frm
            sh_exp = randsample(size(exp_cond, 2), size(exp_cond, 2))'; %shuffle trials for bootstrapping
            exp_ind = sh_exp(ind);
            nexp_cond(timepoint, :) = exp_cond(timepoint, exp_ind(1, :)); %resampled expected trials
            nexp_avg(timepoint, :) = nanmean(nexp_cond(timepoint, :), 2);
            nun_cond(timepoint, :) = un_cond(timepoint, ind(1, :)); %resampled unexpected trials
            nun_avg(timepoint, :) = nanmean(nun_cond(timepoint, :), 2);
            %ndiff_expun_avg(timepoint, :) = nun_avg(timepoint, :) - nexp_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
            %of resampled expected and unexpected trials at each timepoint
            ndiff_expun_avg(timepoint, :) = nexp_avg(timepoint, :) - nun_avg(timepoint, :);
        end
        nexp_cond_b1{nboot} = nexp_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        nun_cond_b1{nboot} = nun_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nexp_avg_b1(:, nboot) = nexp_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        nun_avg_b1(:, nboot) = nun_avg;
        ndiff_expun_avg_b1(:, nboot) = ndiff_expun_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    
    %computing CI
    for tp = 1:frm
        resp.ndiff_expun_b1_CI(:, tp) = prctile(ndiff_expun_avg_b1(tp, :), [2.5 97.5]);
        resp.nexp_b1_CI(:, tp) = prctile(nexp_avg_b1(tp, :), [2.5 97.5]);
        resp.nun_b1_CI(:, tp) = prctile(nun_avg_b1(tp, :), [2.5 97.5]);
        resp.p_exp_b1(:, tp) = min((2*min(numel(find(ndiff_expun_avg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(ndiff_expun_avg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(ndiff_expun_avg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(ndiff_expun_avg_b1(tp, :) > 0))./numboot))));
    end
    
    %%%%%%%%%%%%%%%%%%%resp_attention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nboot = 1:numboot
        clear nfoc_cond ndiv_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
        
        foc_cond = data.resp_foc1{ss}; %all foc trials in b1 of that subj
        div_cond = data.resp_div1{ss}; %all div trials in b1 of that subj
        
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
        nfoc_cond_b1{nboot} = nfoc_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        ndiv_cond_b1{nboot} = ndiv_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nfoc_avg_b1(:, nboot) = nfoc_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        ndiv_avg_b1(:, nboot) = ndiv_avg;
        ndiff_focdiv_avg_b1(:, nboot) = ndiff_focdiv_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
    end
    
    for tp = 1:frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        resp.ndiff_focdiv_b1_CI(:, tp) = prctile(ndiff_focdiv_avg_b1(tp, :), [2.5 97.5]);
        resp.nfoc_b1_CI(:, tp) = prctile(nfoc_avg_b1(tp, :), [2.5 97.5]);
        resp.ndiv_b1_CI(:, tp) = prctile(ndiv_avg_b1(tp, :), [2.5 97.5]);
        resp.p_att_b1(:, tp) = min((2*min(numel(find(ndiff_focdiv_avg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(ndiff_focdiv_avg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(ndiff_focdiv_avg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(ndiff_focdiv_avg_b1(tp, :) > 0))./numboot))));
    end
    
    %%%%%%%%%%%%%%%%%%%resp_coh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nboot = 1:numboot
        clear nhi_cond nlo_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
        
        hi_cond = data.resp_hi1{ss}; %all foc trials in b1 of that subj
        lo_cond = data.resp_lo1{ss}; %all div trials in b1 of that subj
        
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
        nhi_cond_b1{nboot} = nhi_cond; %resampled exp trials; size = frm x size(un_cond, 2)
        nlo_cond_b1{nboot} = nlo_cond; %resampled un trials; size = frm x size(un_cond, 2)
        
        nhi_avg_b1(:, nboot) = nhi_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
        nlo_avg_b1(:, nboot) = nlo_avg;
        ndiff_hilo_avg_b1(:, nboot) = ndiff_hilo_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
        
    end
    
    for tp = 1:frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        resp.ndiff_hilo_b1_CI(:, tp) = prctile(ndiff_hilo_avg_b1(tp, :), [2.5 97.5]);
        resp.nhi_b1_CI(:, tp) = prctile(nhi_avg_b1(tp, :), [2.5 97.5]);
        resp.nlo_b1_CI(:, tp) = prctile(nlo_avg_b1(tp, :), [2.5 97.5]);
        resp.p_coh_b1(:, tp) = min((2*min(numel(find(ndiff_hilo_avg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(ndiff_hilo_avg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(ndiff_hilo_avg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(ndiff_hilo_avg_b1(tp, :) > 0))./numboot))));
    end
    
    
    %%
    %plotting things
    figure(ss); %resp trajectories
    xlabel('Time (ms)')
    ylabel('Response error (deg)')
    %set(n2, 'Fontsize', 13)
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    %     ylim([0 1.0]); YTick([0:0.2:1]);
    %     xlim([0 1500]); XTick([0:500:1500]);
    
    subplot(3, 1, 3) %exp
    %p1 = plot(timex, nanmean(nexp_avg_b1, 2), 'Color', c_exp, 'LineWidth', 3); hold on; %resampled expected trials
    p1 = plot(timex, resp.ndiff_expun_b1_CI(1, :), 'Color', c_pat, 'LineWidth', 1.5); hold on;
    p2 = plot(timex, resp.ndiff_expun_b1_CI(2, :), 'Color', c_pat, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, (-1)*cresp.ndiff_expun_b1_CI(1, :), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p4 = plot(timex, (-1)*cresp.ndiff_expun_b1_CI(2, :), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p5 = plot(timex, (-0.5)* (cresp.ndiff_expun_b1_CI(1, :)+ cresp.ndiff_expun_b1_CI(2, :)), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p3.Color(4) = opac;
    p4.Color(4) = opac;
    title('Resp error: Expected minus unexpected');
    ylabel('Response error (deg)')
    legend([p1 p3], {'Patient', 'Mean controls'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 2) %att
    p1 = plot(timex, resp.ndiff_focdiv_b1_CI(1, :), 'Color', c_pat, 'LineWidth', 1.5); hold on;
    p2 = plot(timex, resp.ndiff_focdiv_b1_CI(2, :), 'Color', c_pat, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, (-1)*cresp.ndiff_focdiv_b1_CI(1, :), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p4 = plot(timex, (-1)*cresp.ndiff_focdiv_b1_CI(2, :), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p5 = plot(timex, (-0.5)* (cresp.ndiff_focdiv_b1_CI(1, :)+ cresp.ndiff_focdiv_b1_CI(2, :)), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p3.Color(4) = opac;
    p4.Color(4) = opac;
    title('Resp error: Foc minus div');
    ylabel('Response error (deg)')
    legend([p1 p3], {'Patient', 'Mean controls'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    subplot(3, 1, 1) %coh
    p1 = plot(timex, resp.ndiff_hilo_b1_CI(1, :), 'Color', c_pat, 'LineWidth', 1.5); hold on;
    p2 = plot(timex, resp.ndiff_hilo_b1_CI(2, :), 'Color', c_pat, 'LineWidth', 1.5); hold on;
    p3 = plot(timex, (-1)*cresp.ndiff_hilo_b1_CI(1, :), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p4 = plot(timex, (-1)*cresp.ndiff_hilo_b1_CI(2, :), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p5 = plot(timex, (-0.5)* (cresp.ndiff_hilo_b1_CI(1, :)+ cresp.ndiff_hilo_b1_CI(2, :)), 'Color', c_con, 'LineWidth', 1.5); hold on;
    p3.Color(4) = opac;
    p4.Color(4) = opac;
    title('Resp error: High minus low coh');
    ylabel('Response error (deg)')
    legend([p1 p3], {'Patient', 'Mean controls'});
    set(gca, 'FontSize', 11, 'fontWeight', 'bold')
    
    %setting subplot properties
    AX = subplot(3, 1, 1);
    AY = subplot(3, 1, 2);
    AZ = subplot(3, 1, 3);
    set([AX AY AZ], 'YLim', [-35 30]);
    set([AX AY AZ], 'YTick', [-35:10:30]);
    set([AX AY AZ], 'XLim', [0 1500]);
    set([AX AY AZ], 'XTick', [0:500:1500]);
end
toc

