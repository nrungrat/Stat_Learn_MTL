% This script compares response onsets between patients and controls
close all; clear all;

cc = load('controls_data.mat')
pp = load('patients_data.mat')

numboot = 1000;
numboot2 = 100;

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
noise = 0.008;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Resp Trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1A) resp trajectories as a function of expectation
cc.traj.onset_exp = nan(1, numboot2);
pp.traj.onset_exp = nan(1, numboot2);
for nboot2 = 1:numboot2 %numboot
    nboot2
    for nboot = 1:numboot
        %controls
        fprintf('controls: resp traj vs exp_bin1: %s\n', num2str(nboot))
        sub_ind = randperm(numel(cc.allsub), numel(pp.allsub));
        %sub_ind = 1:numel(cc.allsub);
        %sub_ind = randi(numel(cc.allsub),[1, numel(cc.allsub)]); %resampling across subj
        for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
            clear nexp_cond nun_cond ndiff_expun_avg %doing this b/c #of trials in each bin vary across subj
            
            exp_cond = cc.data.dist_exp1{sub_ind(ss)}; %all exp trials in b1 of that subj
            un_cond = cc.data.dist_un1{sub_ind(ss)}; %all un trials in b1 of that subj
            
            %%% 2 ways of doing this i think: 1) fixed # across subj determined by the
            %%% min # of un trials across subj; 2) the number of un of each subj
            %%% trying 2) first
            ind = randi(size(un_cond, 2), [1, size(un_cond, 2)]); % this is an index matrix to sample
            % out XX trials (with replacement) each from expected and
            % unexpected condition
            for timepoint = 1:cc.frm
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
            %         %onset
            %         dif = (ndiff_expun_avg(2:end))-(ndiff_expun_avg(1:end-1));
            %         indd = find(abs(dif) > noise, 1);
            %
            %         if isempty(indd)
            %             cc.onset_expun(:, ss) = nan;
            %         else
            %             cc.onset_expun(:, ss) = cc.timex(indd); %the 'first' point that traj deviates from baseline
            %         end
        end
        
        cc.traj.nexp_savg_b1(:, nboot) = nanmean(nexp_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        cc.traj.nun_savg_b1(:, nboot) = nanmean(nun_avg_b1, 2);
        cc.traj.ndiff_expun_savg_b1(:, nboot) = nanmean(ndiff_expun_avg_b1, 2);
        
        %     %find onset
        %     onset_ind = cc.timex(find(abs(nanmean(ndiff_expun_avg_b1, 2)) > noise));
        %     cc.traj.onset_expun(:, nboot) = onset_ind(find(onset_ind > 500, 1));
        %
        %%----------------------patients---------------------------------------
        fprintf('patients: resp traj vs exp_bin1: %s\n', num2str(nboot))
        sub_ind = 1:numel(pp.allsub);
        %sub_ind = randi(numel(pp.allsub),[1, numel(pp.allsub)]); %resampling across subj
        for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
            clear nexp_cond nun_cond ndiff_expun_avg %doing this b/c #of trials in each bin vary across subj
            
            exp_cond = pp.data.dist_exp1{sub_ind(ss)}; %all exp trials in b1 of that subj
            un_cond = pp.data.dist_un1{sub_ind(ss)}; %all un trials in b1 of that subj
            
            %%% 2 ways of doing this i think: 1) fixed # across subj determined by the
            %%% min # of un trials across subj; 2) the number of un of each subj
            %%% trying 2) first
            ind = randi(size(un_cond, 2), [1, size(un_cond, 2)]); % this is an index matrix to sample
            % out XX trials (with replacement) each from expected and
            % unexpected condition
            for timepoint = 1:pp.frm
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
            %         %onset
            %         dif = (ndiff_expun_avg(2:end))-(ndiff_expun_avg(1:end-1));
            %         indd = find(abs(dif) > noise, 1);
            %
            %         if isempty(indd)
            %             pp.onset_expun(:, ss) = nan;
            %         else
            %             pp.onset_expun(:, ss) = pp.timex(indd); %the 'first' point that traj deviates from baseline
            %         end
        end
        
        pp.traj.nexp_savg_b1(:, nboot) = nanmean(nexp_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        pp.traj.nun_savg_b1(:, nboot) = nanmean(nun_avg_b1, 2);
        pp.traj.ndiff_expun_savg_b1(:, nboot) = nanmean(ndiff_expun_avg_b1, 2);
        
        %     %%find onset
        %     onset_ind = pp.timex(find(abs(nanmean(ndiff_expun_avg_b1, 2)) > noise));
        %     pp.traj.onset_expun(:, nboot) = onset_ind(find(onset_ind > 500, 1));
    end
    
    %diff: computing significance
    for tp = 1:cc.frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        cc.traj.ndiff_expun_b1_CI(:, tp) = prctile(cc.traj.ndiff_expun_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.nexp_b1_CI(:, tp) = prctile(cc.traj.nexp_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.nun_b1_CI(:, tp) = prctile(cc.traj.nun_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.p_exp_b1(:, tp) = min((2*min(numel(find(cc.traj.ndiff_expun_savg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(cc.traj.ndiff_expun_savg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(cc.traj.ndiff_expun_savg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(cc.traj.ndiff_expun_savg_b1(tp, :) > 0))./numboot))));
    end
    for tp = 1:pp.frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        pp.traj.ndiff_expun_b1_CI(:, tp) = prctile(pp.traj.ndiff_expun_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.nexp_b1_CI(:, tp) = prctile(pp.traj.nexp_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.nun_b1_CI(:, tp) = prctile(pp.traj.nun_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.p_exp_b1(:, tp) = min((2*min(numel(find(pp.traj.ndiff_expun_savg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(pp.traj.ndiff_expun_savg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(pp.traj.ndiff_expun_savg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(pp.traj.ndiff_expun_savg_b1(tp, :) > 0))./numboot))));
    end
    
    %
    % %onset con vs patients: computing sig
    %     traj.onset_expun_con_vs_pat = pp.traj.onset_expun - cc.traj.onset_expun;
    %     traj.onset_expun_b1_CI = prctile(traj.onset_expun_con_vs_pat, [2.5 97.5]);
    % %     traj.nexp_b1_CI(:, tp) = prctile(traj.nexp_savg_b1(tp, :), [2.5 97.5]);
    % %     traj.nun_b1_CI(:, tp) = prctile(traj.nun_savg_b1(tp, :), [2.5 97.5]);
    %     traj.onset_exp_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_expun_con_vs_pat > 0))./numboot, ...
    %         1-(numel(find(traj.onset_expun_con_vs_pat < 0))./numboot))), ...
    %     (2*min(numel(find(traj.onset_expun_con_vs_pat < 0))./numboot, ...
    %         1-(numel(find(traj.onset_expun_con_vs_pat > 0))./numboot))));
    %
    
    %find where the sig effects of exp starts
    cc_onset_ind = find(cc.traj.p_exp_b1< 0.05, 1);
    pp_onset_ind = find(pp.traj.p_exp_b1< 0.05, 1);
    
    if ~isempty(cc_onset_ind)
        cc.traj.onset_exp(nboot2) =  cc.timex(cc_onset_ind);
    end
    if ~isempty(pp_onset_ind)
        pp.traj.onset_exp(nboot2) =  pp.timex(pp_onset_ind);
    end
    clear cc_onset_ind pp_onset_ind
end

%onset con vs patients: computing sig
con_onset_exp = cc.traj.onset_exp;
pat_onset_exp = pp.traj.onset_exp;
traj.onset_expun_con_vs_pat = pat_onset_exp - con_onset_exp;
traj.onset_expun_b1_CI = prctile(traj.onset_expun_con_vs_pat, [2.5 97.5]);
traj.onset_exp_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_expun_con_vs_pat > 0))./numboot2, ...
    1-(numel(find(traj.onset_expun_con_vs_pat < 0))./numboot2))), ...
    (2*min(numel(find(traj.onset_expun_con_vs_pat < 0))./numboot2, ...
    1-(numel(find(traj.onset_expun_con_vs_pat > 0))./numboot2))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2A) resp trajectories as a function of attention
cc.traj.onset_att = nan(1, numboot2);
pp.traj.onset_att = nan(1, numboot2);
for nboot2 = 1:numboot2 %numboot
    nboot2
    for nboot = 1:numboot
        %fprintf('controls: resp traj vs att_bin1: %s\n', num2str(nboot))
        sub_ind = randperm(numel(cc.allsub), numel(pp.allsub));
        %sub_ind = 1:numel(cc.allsub);
        %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
        for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
            clear nfoc_cond ndiv_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
            
            foc_cond = cc.data.dist_foc1{sub_ind(ss)}; %all foc trials in b1 of that subj
            div_cond = cc.data.dist_div1{sub_ind(ss)}; %all div trials in b1 of that subj
            
            %use the smaller cond (foc or div) to find index
            min_focdiv = min(size(foc_cond, 2), size(div_cond, 2));
            ind = randi(min_focdiv, [1, min_focdiv]); % this is an index matrix to sample
            
            for timepoint = 1:cc.frm
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
        %traj.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        %traj.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
        %traj.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
        
        cc.traj.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        cc.traj.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
        cc.traj.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
        
        %     %find onset
        %     onset_ind = cc.timex(find(abs(nanmean(ndiff_focdiv_avg_b1, 2)) > noise));
        %     cc.traj.onset_focdiv(:, nboot) = onset_ind(find(onset_ind > 500, 1));
        %
        % --------------------------- patients ------------------------------------
        
        %fprintf('patients: resp traj vs att_bin1: %s\n', num2str(nboot))
        sub_ind = 1:numel(pp.allsub);
        %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
        for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
            clear nfoc_cond ndiv_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
            
            foc_cond = pp.data.dist_foc1{sub_ind(ss)}; %all foc trials in b1 of that subj
            div_cond = pp.data.dist_div1{sub_ind(ss)}; %all div trials in b1 of that subj
            
            %use the smaller cond (foc or div) to find index
            min_focdiv = min(size(foc_cond, 2), size(div_cond, 2));
            ind = randi(min_focdiv, [1, min_focdiv]); % this is an index matrix to sample
            
            for timepoint = 1:pp.frm
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
        %traj.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        %traj.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
        %traj.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
        
        pp.traj.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        pp.traj.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
        pp.traj.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
        
        %%find onset
        %onset_ind = pp.timex(find(abs(nanmean(ndiff_focdiv_avg_b1, 2)) > noise));
        %pp.traj.onset_focdiv(:, nboot) = onset_ind(find(onset_ind > 500, 1));
    end
    
    % for tp = 1:cc.frm %loop over timepoints to find average across iterations
    %     %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    %     traj.ndiff_focdiv_b1_CI(:, tp) = prctile(traj.ndiff_focdiv_savg_b1(tp, :), [2.5 97.5]);
    %     traj.nfoc_b1_CI(:, tp) = prctile(traj.nfoc_savg_b1(tp, :), [2.5 97.5]);
    %     traj.ndiv_b1_CI(:, tp) = prctile(traj.ndiv_savg_b1(tp, :), [2.5 97.5]);
    %     traj.p_att_b1(:, tp) = min((2*min(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot, ...
    %         1-(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot))), ...
    %     (2*min(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot, ...
    %         1-(numel(find(traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot))));
    % end
    
    % %onset con vs patients: computing sig
    % traj.onset_focdiv_con_vs_pat = pp.traj.onset_focdiv - cc.traj.onset_focdiv;
    % traj.onset_focdiv_b1_CI = prctile(traj.onset_focdiv_con_vs_pat, [2.5 97.5]);
    % %     traj.nexp_b1_CI(:, tp) = prctile(traj.nexp_savg_b1(tp, :), [2.5 97.5]);
    % %     traj.nun_b1_CI(:, tp) = prctile(traj.nun_savg_b1(tp, :), [2.5 97.5]);
    % traj.onset_att_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_focdiv_con_vs_pat > 0))./numboot, ...
    %     1-(numel(find(traj.onset_focdiv_con_vs_pat < 0))./numboot))), ...
    %     (2*min(numel(find(traj.onset_focdiv_con_vs_pat < 0))./numboot, ...
    %     1-(numel(find(traj.onset_focdiv_con_vs_pat > 0))./numboot))));
    
    %diff: computing significance
    for tp = 1:cc.frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        cc.traj.ndiff_focdiv_b1_CI(:, tp) = prctile(cc.traj.ndiff_focdiv_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.nfoc_b1_CI(:, tp) = prctile(cc.traj.nfoc_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.ndiv_b1_CI(:, tp) = prctile(cc.traj.ndiv_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.p_att_b1(:, tp) = min((2*min(numel(find(cc.traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(cc.traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(cc.traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(cc.traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot))));
    end
    for tp = 1:pp.frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        pp.traj.ndiff_focdiv_b1_CI(:, tp) = prctile(pp.traj.ndiff_focdiv_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.nfoc_b1_CI(:, tp) = prctile(pp.traj.nfoc_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.ndiv_b1_CI(:, tp) = prctile(pp.traj.ndiv_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.p_att_b1(:, tp) = min((2*min(numel(find(pp.traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(pp.traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(pp.traj.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(pp.traj.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot))));
    end
    
    %find where the sig effects of exp starts
    cc_onset_ind = find(cc.traj.p_att_b1< 0.05, 1);
    pp_onset_ind = find(pp.traj.p_att_b1< 0.05, 1);
    
    if ~isempty(cc_onset_ind)
        cc.traj.onset_att(nboot2) =  cc.timex(cc_onset_ind);
    end
    if ~isempty(pp_onset_ind)
        pp.traj.onset_att(nboot2) =  pp.timex(pp_onset_ind);
    end
    clear cc_onset_ind pp_onset_ind
end

% %onset con vs patients: computing sig
% traj.onset_focdiv_con_vs_pat = pp.traj.onset_focdiv - cc.traj.onset_focdiv;
% traj.onset_focdiv_b1_CI = prctile(traj.onset_focdiv_con_vs_pat, [2.5 97.5]);
% %     traj.nexp_b1_CI(:, tp) = prctile(traj.nexp_savg_b1(tp, :), [2.5 97.5]);
% %     traj.nun_b1_CI(:, tp) = prctile(traj.nun_savg_b1(tp, :), [2.5 97.5]);
% traj.onset_att_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_focdiv_con_vs_pat > 0))./numboot, ...
%     1-(numel(find(traj.onset_focdiv_con_vs_pat < 0))./numboot))), ...
%     (2*min(numel(find(traj.onset_focdiv_con_vs_pat < 0))./numboot, ...
%     1-(numel(find(traj.onset_focdiv_con_vs_pat > 0))./numboot))));

%onset con vs patients: computing sig
con_onset_att = cc.traj.onset_att;
pat_onset_att = pp.traj.onset_att;
traj.onset_focdiv_con_vs_pat = pat_onset_att - con_onset_att;
traj.onset_focdiv_b1_CI = prctile(traj.onset_focdiv_con_vs_pat, [2.5 97.5]);
traj.onset_att_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_focdiv_con_vs_pat > 0))./numboot2, ...
    1-(numel(find(traj.onset_focdiv_con_vs_pat < 0))./numboot2))), ...
    (2*min(numel(find(traj.onset_focdiv_con_vs_pat < 0))./numboot2, ...
    1-(numel(find(traj.onset_focdiv_con_vs_pat > 0))./numboot2))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3A) resp trajectories as a function of coherence
cc.traj.onset_coh = nan(1, numboot2);
pp.traj.onset_coh = nan(1, numboot2);
for nboot2 = 1:numboot2 %numboot
    nboot2
    for nboot = 1:numboot
        %fprintf('controls: resp traj vs coh_bin1: %s\n', num2str(nboot))
        sub_ind = randperm(numel(cc.allsub), numel(pp.allsub));
        %sub_ind = 1:numel(cc.allsub);
        %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
        for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
            clear nhi_cond nlo_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
            
            hi_cond = cc.data.dist_hi1{sub_ind(ss)}; %all foc trials in b1 of that subj
            lo_cond = cc.data.dist_lo1{sub_ind(ss)}; %all div trials in b1 of that subj
            
            %use the smaller cond (foc or div) to find index
            min_hilo = min(size(hi_cond, 2), size(lo_cond, 2));
            ind = randi(min_hilo, [1, min_hilo]); % this is an index matrix to sample
            
            for timepoint = 1:cc.frm
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
        %     traj.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        %     traj.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
        %     traj.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
        
        cc.traj.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        cc.traj.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
        cc.traj.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
        
        %     %find onset
        %     onset_ind = cc.timex(find(abs(nanmean(ndiff_hilo_avg_b1, 2)) > noise));
        %     cc.traj.onset_hilo(:, nboot) = onset_ind(find(onset_ind > 500, 1));
        %
        %---------------------------- patients-------------------------------
        %fprintf('patients: resp traj vs coh_bin1: %s\n', num2str(nboot))
        sub_ind = 1:numel(pp.allsub);
        %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
        for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
            clear nhi_cond nlo_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
            
            hi_cond = pp.data.dist_hi1{sub_ind(ss)}; %all foc trials in b1 of that subj
            lo_cond = pp.data.dist_lo1{sub_ind(ss)}; %all div trials in b1 of that subj
            
            %use the smaller cond (foc or div) to find index
            min_hilo = min(size(hi_cond, 2), size(lo_cond, 2));
            ind = randi(min_hilo, [1, min_hilo]); % this is an index matrix to sample
            
            for timepoint = 1:pp.frm
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
        %     traj.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        %     traj.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
        %     traj.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
        
        pp.traj.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
        pp.traj.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
        pp.traj.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
        
        %     %find onset
        %     onset_ind = pp.timex(find(abs(nanmean(ndiff_hilo_avg_b1, 2)) > noise));
        %     pp.traj.onset_hilo(:, nboot) = onset_ind(find(onset_ind > 500, 1));
        %
    end
    
    % for tp = 1:cc.frm %loop over timepoints to find average across iterations
    %     %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
    %     traj.ndiff_hilo_b1_CI(:, tp) = prctile(traj.ndiff_hilo_savg_b1(tp, :), [2.5 97.5]);
    %     traj.nhi_b1_CI(:, tp) = prctile(traj.nhi_savg_b1(tp, :), [2.5 97.5]);
    %     traj.nlo_b1_CI(:, tp) = prctile(traj.nlo_savg_b1(tp, :), [2.5 97.5]);
    %     traj.p_coh_b1(:, tp) = min((2*min(numel(find(traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot, ...
    %         1-(numel(find(traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot))), ...
    %     (2*min(numel(find(traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot, ...
    %         1-(numel(find(traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot))));
    % end
    
    % %onset con vs patients: computing sig
    % traj.onset_hilo_con_vs_pat = pp.traj.onset_hilo - cc.traj.onset_hilo;
    % traj.onset_hilo_b1_CI = prctile(traj.onset_hilo_con_vs_pat, [2.5 97.5]);
    % %     traj.nexp_b1_CI(:, tp) = prctile(traj.nexp_savg_b1(tp, :), [2.5 97.5]);
    % %     traj.nun_b1_CI(:, tp) = prctile(traj.nun_savg_b1(tp, :), [2.5 97.5]);
    % traj.onset_coh_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_hilo_con_vs_pat > 0))./numboot, ...
    %     1-(numel(find(traj.onset_hilo_con_vs_pat < 0))./numboot))), ...
    %     (2*min(numel(find(traj.onset_hilo_con_vs_pat < 0))./numboot, ...
    %     1-(numel(find(traj.onset_hilo_con_vs_pat > 0))./numboot))));
    
    %diff: computing significance
    for tp = 1:cc.frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        cc.traj.ndiff_hilo_b1_CI(:, tp) = prctile(cc.traj.ndiff_hilo_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.nhi_b1_CI(:, tp) = prctile(cc.traj.nhi_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.nlo_b1_CI(:, tp) = prctile(cc.traj.nlo_savg_b1(tp, :), [2.5 97.5]);
        cc.traj.p_coh_b1(:, tp) = min((2*min(numel(find(cc.traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(cc.traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(cc.traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(cc.traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot))));
    end
    for tp = 1:pp.frm %loop over timepoints to find average across iterations
        %display([nanmean(traj.ndiff_savg_b1(tp, :)), prctile(traj.ndiff_savg_b1(tp, :), [2.5 97.5])]);
        pp.traj.ndiff_hilo_b1_CI(:, tp) = prctile(pp.traj.ndiff_hilo_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.nhi_b1_CI(:, tp) = prctile(pp.traj.nhi_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.nlo_b1_CI(:, tp) = prctile(pp.traj.nlo_savg_b1(tp, :), [2.5 97.5]);
        pp.traj.p_coh_b1(:, tp) = min((2*min(numel(find(pp.traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot, ...
            1-(numel(find(pp.traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot))), ...
            (2*min(numel(find(pp.traj.ndiff_hilo_savg_b1(tp, :) < 0))./numboot, ...
            1-(numel(find(pp.traj.ndiff_hilo_savg_b1(tp, :) > 0))./numboot))));
    end
    
    %find where the sig effects of exp starts
    cc_onset_ind = find(cc.traj.p_coh_b1< 0.05, 1);
    pp_onset_ind = find(pp.traj.p_coh_b1< 0.05, 1);
    
    if ~isempty(cc_onset_ind)
        cc.traj.onset_coh(nboot2) =  cc.timex(cc_onset_ind);
    end
    if ~isempty(pp_onset_ind)
        pp.traj.onset_coh(nboot2) =  pp.timex(pp_onset_ind);
    end
    clear cc_onset_ind pp_onset_ind
end
%
con_onset_coh = cc.traj.onset_coh;
pat_onset_coh = pp.traj.onset_coh;
traj.onset_hilo_con_vs_pat = pat_onset_coh - con_onset_coh;
traj.onset_hilo_b1_CI = prctile(traj.onset_hilo_con_vs_pat, [2.5 97.5]);
traj.onset_coh_con_vs_pat_p_b1 = min((2*min(numel(find(traj.onset_hilo_con_vs_pat > 0))./numboot2, ...
    1-(numel(find(traj.onset_hilo_con_vs_pat < 0))./numboot2))), ...
    (2*min(numel(find(traj.onset_hilo_con_vs_pat < 0))./numboot2, ...
    1-(numel(find(traj.onset_hilo_con_vs_pat > 0))./numboot2))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% %%%%%%%%%%%%%%% Resp Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 4A) resp error as a function of expectation
% % bin 1
% for nboot = 1:numboot
%     fprintf('resp error vs exp_bin1: %s\n', num2str(nboot))
%     sub_ind = 1:numel(allsub);
%     %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
%     for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
%         clear nexp_cond nun_cond ndiff_expun_avg %doing this b/c #of trials in each bin vary across subj
%
%         exp_cond = data.resp_exp1{sub_ind(ss)}; %all exp trials in b1 of that subj
%         un_cond = data.resp_un1{sub_ind(ss)}; %all un trials in b1 of that subj
%
%         %%% 2 ways of doing this i think: 1) fixed # across subj determined by the
%         %%% min # of un trials across subj; 2) the number of un of each subj
%         %%% trying 2) first
%         ind = randi(size(un_cond, 2), [1, size(un_cond, 2)]); % this is an index matrix to sample
%         % out XX trials (with replacement) each from expected and
%         % unexpected condition
%         for timepoint = 1:cc.frm
%             sh_exp = randsample(size(exp_cond, 2), size(exp_cond, 2))'; %shuffle trials for bootstrapping
%             exp_ind = sh_exp(ind);
%             nexp_cond(timepoint, :) = exp_cond(timepoint, exp_ind(1, :)); %resampled expected trials
%             nexp_avg(timepoint, :) = nanmean(nexp_cond(timepoint, :), 2);
%             nun_cond(timepoint, :) = un_cond(timepoint, ind(1, :)); %resampled unexpected trials
%             nun_avg(timepoint, :) = nanmean(nun_cond(timepoint, :), 2);
%             ndiff_expun_avg(timepoint, :) = nun_avg(timepoint, :) - nexp_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
%             %of resampled expected and unexpected trials at each timepoint
%         end
%         nexp_cond_b1{ss} = nexp_cond; %resampled exp trials; size = frm x size(un_cond, 2)
%         nun_cond_b1{ss} = nun_cond; %resampled un trials; size = frm x size(un_cond, 2)
%
%         nexp_avg_b1(:, ss) = nexp_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
%         nun_avg_b1(:, ss) = nun_avg;
%         ndiff_expun_avg_b1(:, ss) = ndiff_expun_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
%     end
%     resp.nexp_savg_b1(:, nboot) = nanmean(nexp_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
%     resp.nun_savg_b1(:, nboot) = nanmean(nun_avg_b1, 2);
%     resp.ndiff_expun_savg_b1(:, nboot) = nanmean(ndiff_expun_avg_b1, 2);
% end
%
% for tp = 1:frm %loop over timepoints to find average across iterations
%     %display([nanmean(resp.ndiff_savg_b1(tp, :)), prctile(resp.ndiff_savg_b1(tp, :), [2.5 97.5])]);
%     resp.ndiff_expun_b1_CI(:, tp) = prctile(resp.ndiff_expun_savg_b1(tp, :), [2.5 97.5]);
%     resp.nexp_b1_CI(:, tp) = prctile(resp.nexp_savg_b1(tp, :), [2.5 97.5]);
%     resp.nun_b1_CI(:, tp) = prctile(resp.nun_savg_b1(tp, :), [2.5 97.5]);
%     resp.p_exp_b1(:, tp) = min((2*min(numel(find(resp.ndiff_expun_savg_b1(tp, :) > 0))./numboot, ...
%         1-(numel(find(resp.ndiff_expun_savg_b1(tp, :) < 0))./numboot))), ...
%     (2*min(numel(find(resp.ndiff_expun_savg_b1(tp, :) < 0))./numboot, ...
%         1-(numel(find(resp.ndiff_expun_savg_b1(tp, :) > 0))./numboot))));
%
%
% end
%
% %% 5A) resp error as a function of attention
% % bin 1
% for nboot = 1:numboot
%     fprintf('resp error vs att_bin1: %s\n', num2str(nboot))
%     sub_ind = 1:numel(allsub);
%     %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
%     for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
%         clear nfoc_cond ndiv_cond ndiff_focdiv_avg %doing this b/c #of trials in each bin vary across subj
%
%         foc_cond = data.resp_foc1{sub_ind(ss)}; %all foc trials in b1 of that subj
%         div_cond = data.resp_div1{sub_ind(ss)}; %all div trials in b1 of that subj
%
%         %use the smaller cond (foc or div) to find index
%         min_focdiv = min(size(foc_cond, 2), size(div_cond, 2));
%         ind = randi(min_focdiv, [1, min_focdiv]); % this is an index matrix to sample
%
%         for timepoint = 1:frm
%             if size(foc_cond, 2) < size(div_cond, 2)
%                sh_div = randsample(size(div_cond, 2), size(div_cond, 2))'; %shuffle trials for bootstrapping
%                div_ind = sh_div(ind);
%                  nfoc_cond(timepoint, :) = foc_cond(timepoint, ind(1, :)); %resampled expected trials
%             ndiv_cond(timepoint, :) = div_cond(timepoint, div_ind(1, :)); %resampled unexpected trials
%
%             else
%             sh_foc = randsample(size(foc_cond, 2), size(foc_cond, 2))'; %shuffle trials for bootstrapping
%             foc_ind = sh_foc(ind);
%             nfoc_cond(timepoint, :) = foc_cond(timepoint, foc_ind(1, :)); %resampled expected trials
%             ndiv_cond(timepoint, :) = div_cond(timepoint, ind(1, :)); %resampled unexpected trials
%             end
%             nfoc_avg(timepoint, :) = nanmean(nfoc_cond(timepoint, :), 2);
%             ndiv_avg(timepoint, :) = nanmean(ndiv_cond(timepoint, :), 2);
%             ndiff_focdiv_avg(timepoint, :) = ndiv_avg(timepoint, :) - nfoc_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
%             %of resampled expected and unexpected trials at each timepoint
%         end
%         nfoc_cond_b1{ss} = nfoc_cond; %resampled exp trials; size = frm x size(un_cond, 2)
%         ndiv_cond_b1{ss} = ndiv_cond; %resampled un trials; size = frm x size(un_cond, 2)
%
%         nfoc_avg_b1(:, ss) = nfoc_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
%         ndiv_avg_b1(:, ss) = ndiv_avg;
%         ndiff_focdiv_avg_b1(:, ss) = ndiff_focdiv_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
%     end
%     resp.nfoc_savg_b1(:, nboot) = nanmean(nfoc_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
%     resp.ndiv_savg_b1(:, nboot) = nanmean(ndiv_avg_b1, 2);
%     resp.ndiff_focdiv_savg_b1(:, nboot) = nanmean(ndiff_focdiv_avg_b1, 2);
% end
%
% for tp = 1:frm %loop over timepoints to find average across iterations
%     %display([nanmean(resp.ndiff_savg_b1(tp, :)), prctile(resp.ndiff_savg_b1(tp, :), [2.5 97.5])]);
%     resp.ndiff_focdiv_b1_CI(:, tp) = prctile(resp.ndiff_focdiv_savg_b1(tp, :), [2.5 97.5]);
%     resp.nfoc_b1_CI(:, tp) = prctile(resp.nfoc_savg_b1(tp, :), [2.5 97.5]);
%     resp.ndiv_b1_CI(:, tp) = prctile(resp.ndiv_savg_b1(tp, :), [2.5 97.5]);
%     resp.p_att_b1(:, tp) = min((2*min(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot, ...
%         1-(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot))), ...
%     (2*min(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) < 0))./numboot, ...
%         1-(numel(find(resp.ndiff_focdiv_savg_b1(tp, :) > 0))./numboot))));
%
% end
%
% %% 6A) resp trajectories as a function of coherence
% % bin 1
% for nboot = 1:numboot
%     fprintf('resp error vs coh_bin1: %s\n', num2str(nboot))
%     sub_ind = 1:numel(allsub);
%     %sub_ind = randi(numel(allsub),[1, numel(allsub)]); %resampling across subj
%     for ss = 1:numel(sub_ind) %go through subject's data based on subject index (sub_ind)
%         clear nhi_cond nlo_cond ndiff_hilo_avg %doing this b/c #of trials in each bin vary across subj
%
%         hi_cond = data.resp_hi1{sub_ind(ss)}; %all foc trials in b1 of that subj
%         lo_cond = data.resp_lo1{sub_ind(ss)}; %all div trials in b1 of that subj
%
%         %use the smaller cond (foc or div) to find index
%         min_hilo = min(size(hi_cond, 2), size(lo_cond, 2));
%         ind = randi(min_hilo, [1, min_hilo]); % this is an index matrix to sample
%
%         for timepoint = 1:frm
%             if size(hi_cond, 2) < size(lo_cond, 2)
%                sh_lo = randsample(size(lo_cond, 2), size(lo_cond, 2))'; %shuffle trials for bootstrapping
%                lo_ind = sh_lo(ind);
%                  nhi_cond(timepoint, :) = hi_cond(timepoint, ind(1, :)); %resampled expected trials
%             nlo_cond(timepoint, :) = lo_cond(timepoint, lo_ind(1, :)); %resampled unexpected trials
%
%             else
%             sh_hi = randsample(size(hi_cond, 2), size(hi_cond, 2))'; %shuffle trials for bootstrapping
%             hi_ind = sh_hi(ind);
%             nhi_cond(timepoint, :) = hi_cond(timepoint, hi_ind(1, :)); %resampled expected trials
%             nlo_cond(timepoint, :) = lo_cond(timepoint, ind(1, :)); %resampled unexpected trials
%             end
%             nhi_avg(timepoint, :) = nanmean(nhi_cond(timepoint, :), 2);
%             nlo_avg(timepoint, :) = nanmean(nlo_cond(timepoint, :), 2);
%             %ndiff_hilo_avg(timepoint, :) = nhi_avg(timepoint, :) - nlo_avg(timepoint, :);%diff scores based on XX (size(ind)) pairs
%             %of resampled expected and unexpected trials at each timepoint
%             ndiff_hilo_avg(timepoint, :) = nlo_avg(timepoint, :) - nhi_avg(timepoint, :);
%         end
%         nhi_cond_b1{ss} = nhi_cond; %resampled exp trials; size = frm x size(un_cond, 2)
%         nlo_cond_b1{ss} = nlo_cond; %resampled un trials; size = frm x size(un_cond, 2)
%
%         nhi_avg_b1(:, ss) = nhi_avg; %avg across all resampled exp trials for each timepoint; size = frm x 1
%         nlo_avg_b1(:, ss) = nlo_avg;
%         ndiff_hilo_avg_b1(:, ss) = ndiff_hilo_avg; %diff between the avg exp and avg un for each timepoint; size = frm x 1
%     end
%     resp.nhi_savg_b1(:, nboot) = nanmean(nhi_avg_b1, 2); %avg across all subjects so we get one value at each timepoint; size = frm x nboot
%     resp.nlo_savg_b1(:, nboot) = nanmean(nlo_avg_b1, 2);
%     resp.ndiff_hilo_savg_b1(:, nboot) = nanmean(ndiff_hilo_avg_b1, 2);
% end
%
% for tp = 1:frm %loop over timepoints to find average across iterations
%     %display([nanmean(resp.ndiff_savg_b1(tp, :)), prctile(resp.ndiff_savg_b1(tp, :), [2.5 97.5])]);
%     resp.ndiff_hilo_b1_CI(:, tp) = prctile(resp.ndiff_hilo_savg_b1(tp, :), [2.5 97.5]);
%     resp.nhi_b1_CI(:, tp) = prctile(resp.nhi_savg_b1(tp, :), [2.5 97.5]);
%     resp.nlo_b1_CI(:, tp) = prctile(resp.nlo_savg_b1(tp, :), [2.5 97.5]);
%     resp.p_coh_b1(:, tp) = min((2*min(numel(find(resp.ndiff_hilo_savg_b1(tp, :) > 0))./numboot, ...
%         1-(numel(find(resp.ndiff_hilo_savg_b1(tp, :) < 0))./numboot))), ...
%     (2*min(numel(find(resp.ndiff_hilo_savg_b1(tp, :) < 0))./numboot, ...
%         1-(numel(find(resp.ndiff_hilo_savg_b1(tp, :) > 0))./numboot))));
% end
%
% %%
% toc

