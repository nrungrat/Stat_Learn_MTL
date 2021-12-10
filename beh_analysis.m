% This script perform behavioral analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
scnt = 0;
%allsub = [11 12 13 14 15 16 17 18 19 21]; %controls
allsub= [30 31 32 33] ; %patients

for s =  1:size(allsub, 2)
    s
    scnt = scnt+1;
    clear pcat
    
    cd (['data/subj' num2str(allsub(s))]);
    
    if allsub(s) < 10
        in_file = dir(['AxeHC_Subj0' num2str(allsub(s)) '*Run*']);
        %size(in_file) = [32x1]; 32 blocks/subj (8 blocks x 4sess)
    else
        in_file = dir(['AxeHC_Subj' num2str(allsub(s)) '*Run*']);
    end
    
    %load the calibration data of each subject for calculation of accuracy, etc.
    if allsub(s) < 10
        load(['AxeHC_Cali_Subj0' num2str(allsub(s)) '.mat.mat']);
    else
        load(['AxeHC_Cali_Subj' num2str(allsub(s)) '.mat.mat']);
    end
    label = {'attcue','tgColor','moco','stimDir','stimDirReal','joyX','joyY','timeeachframe','expectOri'};
    %attcue = attention condition (1--divided, 2--focused)
    %tgColor = target color (1--red, 2--blue)
    %moco = tags for motion coherence (1--0%, 2--low coh, 3--high coh);
    %       size 2*168 for nontg_moco and tg_moco on each trial
    %stimDir = tag assigned to each direction (1-7)
    %stimDirReal = actual direction in rads
    
    %parameters for plotting
    frm = 300;
    fq = 120; %monitor's refresh rate
    timex = 0: 1000/fq : (frm-1)*(1000/fq); %time axis (x)
    c_foc = [0 0.4 0]; %green
    c_div = [0.4 0 0.6]; %purple
    c_exp = [0 0.4 0.8]; %blue
    c_neu = [0.5 0.5 0.5]; %grey
    c_un = [0.8 0 0]; %orange
    c_hi = [0 0 0]; %black
    c_hi = [0.6 0 0.2]; %maroon
    c_lo = [0.4 0.4 0.4]; %grey
    %         diffcolor1 = [0.5 0.5 0.5]; %exp - unexp
    %         diffcolor2 = [0.5 0.5 0.5];
    %         noise = 0.1;
    
    
    for r = 1:size(in_file,1) %looping through all 32 blocks
        r
        clear p
        load(in_file(r).name);
        frm = 300; %# of frames we're interested in
        %target dur = 1500 ms and min ITI = 500 ms so we can set
        %frm = 240 frm (= 2000 ms)
        
        %really only care about the first xx frame
        p.timeeachframe = p.timeeachframe(1:frm, :);
        
        p.joyX = p.joyx(1:frm, :); %x coordinate of the joystick response; [frm(300)xtrials(104)]
        p.joyY = p.joyy(1:frm, :); %y coordinate of the joystick response; [frm(300)xtrials(104)]
        
        for t = 1:size(p.joyX, 2) %looping through all 104 trials
            %distance of joystick movement
            d(t, :) = sqrt(p.joyX(:, t).^2 + p.joyY(:, t).^2);
            %             figure;
            %             plot(d(t, :));
            %find index for max distance
            [~,ind(t)] = max(d(t, :));
            
            %find (x,y) at max distance on that trial
            x_raw(t) = p.joyx(ind(t), t);
            y_raw(t) = p.joyy(ind(t), t);
            
            %correcting for the joystick's square base
            distc = 1; %corrected joystick distance
            %082318
            rampup_ind = 1:ind(t);
            [c index(t)] = min(abs(d(t, rampup_ind) - distc)); %use this instead to make sure it happens before the peak
            
            %[c index(t)] = min(abs(d - distc)); %find the index of of where 'distance' array
            %is closest to 1 (i.e., our
            %corrected distance on the joy
            %base)
            %to see which 'pair'
            %of p.joyx and p.joyy marks the
            %'end' of joy movement i.e.,
            %this is used to later compute
            %for response angle
            
            %find (x,y) at max corrected distance which is ~1 a.u. for each
            %trial
            x(1, t) = p.joyx(index(t), t); %first row doesn't care if the responses happened outside deadline
            y(1, t) = p.joyy(index(t), t);
            
            if index(t) <= frm %response was made before deadline
                %111918
                %if ind(t) <= frm; %try defining late responses in a different way
                x(2, t) = p.joyx(index(t), t); %2nd row says nan for late responses
                y(2, t) = p.joyy(index(t), t);
                
            else %late response
                x(2, t) = NaN;
                y(2, t) = NaN;
            end
        end
        
        %find angle in rad based on (x,y) at max distance
        %size(x) = [2x104], where the 2nd row is nan when response was
        %late
        [ang, disp] = cart2pol(x,y); %[theta, rho] = cart2pol(x,y); disp = displacement = rho
        angles_all(r, :) = wrapToPi(ang(1, :)); %[-pi pi]; doesn't care about trials with late responses
        angles(r, :) = wrapToPi(ang(2, :)); %trials with late responses have ang set to nan
        index_angles(r, :) = index;
        %redundant with angle but just for trying out purpose
        %theta(r, :) = ang;
        %now add rho variable
        %rho(r, :) = disp;
        
        %motion coh of the target on each trial (2--low coh, 3--high coh)
        p.moco = p.moco(p.moco~=1)';
        
        %concatenate
        if r == 1
            for ii = 1:size(label, 2)
                pcat.(label{ii}) = p.(label{ii});
            end
            pcat.angles_all = angles_all(r, :);
            pcat.angles = angles(r, :);
            pcat.index_angles = index_angles(r, :);
            %add rho
            %%%%%%pcat.rho = disp;
            %%concatenate the joystick trajectory angles
            %pcat.angle = angle(r, :);
            
        else
            catindex = [2 2 2 2 2 2 2 2 2];
            for ii = 1:size(label, 2)
                pcat.(label{ii})  = cat(catindex(ii), pcat.(label{ii}), p.(label{ii}));
            end
            %pcat.angle = cat(2, pcat.angle, angle(r, :));
            %pcat.angle = cat(2, pcat.angle, ang); %the final response angle on each trial determined by the end point
            pcat.angles_all = cat(2, pcat.angles_all, angles_all(r, :));
            pcat.angles = cat(2, pcat.angles, angles(r, :));
            pcat.index_angles = cat(2, pcat.index_angles, index_angles(r, :));
            %pcat.rho = cat(2, pcat.rho, disp); %corrected distance
        end
    end
    
    %calibrated angles from calibration sessions; will be later used to
    %compute accuracy (use this instead of the actual presented angles)
    pcat.caliangle = cali.medianangle(pcat.stimDir);
    
    %distance of joystick movement
    pcat.distance0 = sqrt(pcat.joyX.^2 + pcat.joyY.^2);
    
    pcat.goodtrial = ones(1, 1040);
    for tt = 1:size(pcat.distance0, 2) %going through each trial
        tt
        maxd(tt) = max(pcat.distance0(:, tt)); %find max distance for scaling purpose
        %will take this loop below out in a bit; want to just used the
        %distance0 (unscale) then scale only once later by condition etc
        if maxd(tt) > distc %if more than 1 then scale
            pcat.distance(:, tt) = pcat.distance0(:, tt)./maxd(tt);
        else %if less than or equal to 1 already then keep it as is
            pcat.distance(:, tt) = pcat.distance0(:, tt);
        end
        
        %filter out bad trials
        %scenarios: 1) late responses
        %%% throw out data from certain directions that subjects responded inconsistently
        if allsub(s) == 16
            baddir = 5;
            badid = find(pcat.stimDir==baddir); %find index for trials with presented stim is dir 5
            if maxd(tt) < 0.5 || isnan(pcat.angles(:, tt)) || ismember(tt, badid)
                pcat.goodtrial(:, tt) = 0; %index for which trials to be included
            end
        elseif allsub(s) == 12
            baddir = [3 4 5];
            badid = find(ismember(pcat.stimDir, baddir));
            if maxd(tt) < 0.5 || isnan(pcat.angles(:, tt)) || ismember(tt, badid)
                pcat.goodtrial(:, tt) = 0; %index for which trials to be included
            end
        elseif ~ismember(allsub(s), [12, 16])
            if maxd(tt) < 0.5 || isnan(pcat.angles(:, tt))
                pcat.goodtrial(:, tt) = 0; %index for which trials to be included
            end
        end
        
        
        %compute error response on each frame
        [pcat.an0(:, tt), ~] = cart2pol(pcat.joyX(:, tt), pcat.joyY(:, tt)); %response error on each frame
        pcat.an(:, tt) = wrapToPi(pcat.an0(:, tt));
        pcat.angerror(:, tt) =  abs(wrapTo180(rad2deg(pcat.an(:, tt) - repmat(pcat.caliangle(:, tt), [frm, 1]))));
        %pcat.angerror(:, tt) = rad2deg(an(1, tt)) - rad2deg(repmat(pcat.caliangle(1, tt), [1, 1])
        
    end
    
    %% %%%%%%%%%%%%%%%%% response error on each trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %cmap = colormap(parula(7));
    cmap = colormap(parula(6));
    thres = 2; %setting threshold for 'correct' resp to be within 1.5*std
    %pcat.error = wrapTo180(rad2deg(pcat.angles - pcat.caliangle));
    %082318
    pcat.error = abs(wrapTo180(rad2deg(pcat.angles - pcat.caliangle)));
    %alldirs = [0.4488 1.3464 2.2440 3.1416 4.0392 4.9368 5.8344];
    alldirs = [0.8029 2.0595 3.3161 4.5728 5.8294];
    
    %% %%%%%%%%%%%%%% splitting trials into bins based on degree of resp error %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pcat.trlabel = 1:size(pcat.distance, 2); %trial label
    pcat.prior = ones(1, size(pcat.distance, 2))*3; %3 = unexpected
    pcat.prior(find(pcat.stimDir - pcat.expectOri == 0)) = 1; %expected conditions
    pcat.prior(find(isnan(pcat.expectOri) == 1)) = 2; %neutral conditions
    pcat.exp_error = pcat.error(pcat.prior == 1); %performance (response error) on expected trials
    pcat.un_error = pcat.error(pcat.prior == 3); %that on unexpected trials
    pcat.neu_error = pcat.error(pcat.prior == 2); %that on neutral trials
    
    pcat.tr_exp_error = pcat.trlabel(pcat.prior == 1); %trial# for expected trials
    pcat.tr_un_error = pcat.trlabel(pcat.prior == 3); %unexpected trials
    pcat.tr_neu_error = pcat.trlabel(pcat.prior == 2); %neutral trials
    
    %attention effects now
    pcat.foc_error = pcat.error(pcat.attcue == 1);
    pcat.div_error = pcat.error(pcat.attcue == 2);
    
    pcat.tr_foc_error = pcat.trlabel(pcat.attcue == 1); %trial# for focused trials
    pcat.tr_div_error = pcat.trlabel(pcat.attcue == 2); %divided trials
    
    %coherence effects
    pcat.lo_error = pcat.error(pcat.moco == 2);
    pcat.hi_error = pcat.error(pcat.moco == 3);
    
    pcat.tr_lo_error = pcat.trlabel(pcat.moco == 2); %trial# for lo coh trials
    pcat.tr_hi_error = pcat.trlabel(pcat.moco == 3); %high coherent trials
    %
    %     %only use good trials
    %     pcat.distance0_good = pcat.distance0(:, find(pcat.goodtrial==1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now splitting performance on each of the three prior condition trials
    %into 2 bins
    allbin = 1:2;
    opp_thres = 30; %threshold for response in ~oppposite directions in deg
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %expected trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcat.exp_nthbin = nan(1, size(pcat.exp_error, 2)); %an array containing the nth bin label
    exp_opp_ind = find(abs(pcat.exp_error) >= 180-opp_thres & abs(pcat.exp_error) <= 180);
    pcat.exp_nthbin(exp_opp_ind) = allbin(end);
    exp_perf = [(pcat.tr_exp_error)' (abs(pcat.exp_error))']; %using the real trial label & abs resp error
    exp_perf_sorted = sortrows(exp_perf, 2);
    exp_tr_label_sorted = exp_perf_sorted(:, 1); %sorted real trial label
    exp_abs_err_sorted = exp_perf_sorted(:, 2); %sorted absolute error
    
    %find index of resp in the opposite direction range (threshold defined above)
    exp_lastbn_ind = find(exp_abs_err_sorted >= 180-opp_thres & exp_abs_err_sorted <= 180);
    if isempty(exp_lastbn_ind) %if not a single resp was in opposite directions
        exp_bin1_trials = exp_perf_sorted; %keep all trials for binning
    else
        exp_bin1_trials = [exp_tr_label_sorted(1: exp_lastbn_ind-1) exp_abs_err_sorted(1: exp_lastbn_ind-1)]; %remaining trials after taking out 'opposite dir' resp
    end
    exp_binsz = round(size(exp_bin1_trials, 1)./(size(allbin, 2)-1)); %estimate size for the remaining bins (everything but last)
    %exp_splits = [quantile(exp_remain(:, 2), 1/3), quantile(exp_remain(:, 2), 2/3)]; %find the splitting points
    exp_bin2_trials = [exp_tr_label_sorted(exp_lastbn_ind) exp_abs_err_sorted(exp_lastbn_ind)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %unexpected trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcat.un_nthbin = nan(1, size(pcat.un_error, 2));
    un_opp_ind = find(abs(pcat.un_error) >= 180-opp_thres & abs(pcat.un_error) <= 180);
    pcat.un_nthbin(un_opp_ind) = allbin(end);
    un_perf = [(pcat.tr_un_error)' (abs(pcat.un_error))']; %using the real trial label & abs resp error
    un_perf_sorted = sortrows(un_perf, 2);
    un_tr_label_sorted = un_perf_sorted(:, 1); %sorted trial label
    un_abs_err_sorted = un_perf_sorted(:, 2); %sorted absolute error
    
    %find index of resp in the opposite direction range (threshold defined above)
    un_lastbn_ind = find(un_abs_err_sorted >= 180-opp_thres & un_abs_err_sorted <= 180);
    if isempty(un_lastbn_ind) %if not a single resp was in opposite directions
        un_bin1_trials = un_perf_sorted; %keep all trials for binning
    else
        un_bin1_trials = [un_tr_label_sorted(1: un_lastbn_ind-1) un_abs_err_sorted(1: un_lastbn_ind-1)]; %remaining trials after taking out 'opposite dir' resp
    end
    un_binsz = round(size(un_bin1_trials, 1)./(size(allbin, 2)-1)); %estimate size for the remaining bins (everything but last)
    %un_splits = [quantile(un_remain(:, 2), 1/3), quantile(un_remain(:, 2), 2/3)]; %find the splitting points
    un_bin2_trials = [un_tr_label_sorted(un_lastbn_ind) un_abs_err_sorted(un_lastbn_ind)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %focused trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcat.foc_nthbin = nan(1, size(pcat.foc_error, 2)); %an array containing the nth bin label
    foc_opp_ind = find(abs(pcat.foc_error) >= 180-opp_thres & abs(pcat.foc_error) <= 180);
    pcat.foc_nthbin(foc_opp_ind) = allbin(end);
    foc_perf = [(pcat.tr_foc_error)' (abs(pcat.foc_error))']; %using the real trial label & abs resp error
    foc_perf_sorted = sortrows(foc_perf, 2);
    foc_tr_label_sorted = foc_perf_sorted(:, 1); %sorted trial label
    foc_abs_err_sorted = foc_perf_sorted(:, 2); %sorted absolute error
    
    %find index of resp in the opposite direction range (threshold defined above)
    foc_lastbn_ind = find(foc_abs_err_sorted >= 180-opp_thres & foc_abs_err_sorted <= 180);
    if isempty(foc_lastbn_ind) %if not a single resp was in opposite directions
        foc_bin1_trials = foc_perf_sorted; %keep all trials for binning
    else
        foc_bin1_trials = [foc_tr_label_sorted(1: foc_lastbn_ind-1) foc_abs_err_sorted(1: foc_lastbn_ind-1)]; %remaining trials after taking out 'opposite dir' resp
    end
    foc_binsz = round(size(foc_bin1_trials, 1)./(size(allbin, 2)-1)); %estimate size for the remaining bins (everything but last)
    %foc_splits = [quantile(foc_remain(:, 2), 1/3), quantile(foc_remain(:, 2), 2/3)]; %find the splitting points
    foc_bin2_trials = [foc_tr_label_sorted(foc_lastbn_ind) foc_abs_err_sorted(foc_lastbn_ind)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %divided trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcat.div_nthbin = nan(1, size(pcat.div_error, 2)); %an array containing the nth bin label
    div_opp_ind = find(abs(pcat.div_error) >= 180-opp_thres & abs(pcat.div_error) <= 180);
    pcat.div_nthbin(div_opp_ind) = allbin(end);
    div_perf = [(pcat.tr_div_error)' (abs(pcat.div_error))']; %using the real trial label & abs resp error
    div_perf_sorted = sortrows(div_perf, 2);
    div_tr_label_sorted = div_perf_sorted(:, 1); %sorted trial label
    div_abs_err_sorted = div_perf_sorted(:, 2); %sorted absolute error
    
    %find index of resp in the opposite direction range (threshold defined above)
    div_lastbn_ind = find(div_abs_err_sorted >= 180-opp_thres & div_abs_err_sorted <= 180);
    if isempty(div_lastbn_ind) %if not a single resp was in opposite directions
        div_bin1_trials = div_perf_sorted; %keep all trials for binning
    else
        div_bin1_trials = [div_tr_label_sorted(1: div_lastbn_ind-1) div_abs_err_sorted(1: div_lastbn_ind-1)]; %remaining trials after taking out 'opposite dir' resp
    end
    div_binsz = round(size(div_bin1_trials, 1)./(size(allbin, 2)-1)); %estimate size for the remaining bins (everything but last)
    %div_splits = [quantile(div_remain(:, 2), 1/3), quantile(div_remain(:, 2), 2/3)]; %find the splitting points
    div_bin2_trials = [div_tr_label_sorted(div_lastbn_ind) div_abs_err_sorted(div_lastbn_ind)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %low coh trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcat.lo_nthbin = nan(1, size(pcat.lo_error, 2)); %an array containing the nth bin label
    lo_opp_ind = find(abs(pcat.lo_error) >= 180-opp_thres & abs(pcat.lo_error) <= 180);
    pcat.lo_nthbin(lo_opp_ind) = allbin(end);
    lo_perf = [(pcat.tr_lo_error)' (abs(pcat.lo_error))']; %using the real trial label & abs resp error
    lo_perf_sorted = sortrows(lo_perf, 2);
    lo_tr_label_sorted = lo_perf_sorted(:, 1); %sorted trial label
    lo_abs_err_sorted = lo_perf_sorted(:, 2); %sorted absolute error
    
    %find index of resp in the opposite direction range (threshold defined above)
    lo_lastbn_ind = find(lo_abs_err_sorted >= 180-opp_thres & lo_abs_err_sorted <= 180);
    if isempty(lo_lastbn_ind) %if not a single resp was in opposite directions
        lo_bin1_trials = lo_perf_sorted; %keep all trials for binning
    else
        lo_bin1_trials = [lo_tr_label_sorted(1: lo_lastbn_ind-1) lo_abs_err_sorted(1: lo_lastbn_ind-1)]; %remaining trials after taking out 'opposite dir' resp
    end
    lo_binsz = round(size(lo_bin1_trials, 1)./(size(allbin, 2)-1)); %estimate size for the remaining bins (everything but last)
    %lo_splits = [quantile(lo_remain(:, 2), 1/3), quantile(lo_remain(:, 2), 2/3)]; %find the splitting points
    lo_bin2_trials = [lo_tr_label_sorted(lo_lastbn_ind) lo_abs_err_sorted(lo_lastbn_ind)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %hi coh trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcat.hi_nthbin = nan(1, size(pcat.hi_error, 2)); %an array containing the nth bin label
    hi_opp_ind = find(abs(pcat.hi_error) >= 180-opp_thres & abs(pcat.hi_error) <= 180);
    pcat.hi_nthbin(hi_opp_ind) = allbin(end);
    hi_perf = [(pcat.tr_hi_error)' (abs(pcat.hi_error))']; %using the real trial label & abs resp error
    hi_perf_sorted = sortrows(hi_perf, 2);
    hi_tr_label_sorted = hi_perf_sorted(:, 1); %sorted trial label
    hi_abs_err_sorted = hi_perf_sorted(:, 2); %sorted absolute error
    
    %find index of resp in the opposite direction range (threshold defined above)
    hi_lastbn_ind = find(hi_abs_err_sorted >= 180-opp_thres & hi_abs_err_sorted <= 180);
    if isempty(hi_lastbn_ind) %if not a single resp was in opposite directions
        hi_bin1_trials = hi_perf_sorted; %keep all trials for binning
    else
        hi_bin1_trials = [hi_tr_label_sorted(1: hi_lastbn_ind-1) hi_abs_err_sorted(1: hi_lastbn_ind-1)]; %remaining trials after taking out 'opposite dir' resp
    end
    hi_binsz = round(size(hi_bin1_trials, 1)./(size(allbin, 2)-1)); %estimate size for the remaining bins (everything but last)
    %hi_splits = [quantile(hi_remain(:, 2), 1/3), quantile(hi_remain(:, 2), 2/3)]; %find the splitting points
    hi_bin2_trials = [hi_tr_label_sorted(hi_lastbn_ind) hi_abs_err_sorted(hi_lastbn_ind)];
    
    %% filtering out bad trials
    
    badind = find(pcat.goodtrial == 0); %find labels for bad trials
    pcat.distance_good = pcat.distance(:, ~ismember(1:1040, badind));
    pcat.angerror_good = pcat.angerror(:, ~ismember(1:1040, badind));
    all.dist{scnt} = pcat.distance_good;
    all.ang{scnt} = pcat.angerror_good;
    
    %exp/un
    exp_bin1 = exp_bin1_trials(~ismember(exp_bin1_trials(:, 1), badind), :); %only use data from not-bad trials
    exp_bin2 = exp_bin2_trials(~ismember(exp_bin2_trials(:, 1), badind), :);
    un_bin1 = un_bin1_trials(~ismember(un_bin1_trials(:, 1), badind), :); %only use data from not-bad trials
    un_bin2 = un_bin2_trials(~ismember(un_bin2_trials(:, 1), badind), :);
    %foc/div
    foc_bin1 = foc_bin1_trials(~ismember(foc_bin1_trials(:, 1), badind), :); %only use data from not-bad trials
    foc_bin2 = foc_bin2_trials(~ismember(foc_bin2_trials(:, 1), badind), :);
    div_bin1 = div_bin1_trials(~ismember(div_bin1_trials(:, 1), badind), :); %only use data from not-bad trials
    div_bin2 = div_bin2_trials(~ismember(div_bin2_trials(:, 1), badind), :);
    %hi/lo
    hi_bin1 = hi_bin1_trials(~ismember(hi_bin1_trials(:, 1), badind), :); %only use data from not-bad trials
    hi_bin2 = hi_bin2_trials(~ismember(hi_bin2_trials(:, 1), badind), :);
    lo_bin1 = lo_bin1_trials(~ismember(lo_bin1_trials(:, 1), badind), :); %only use data from not-bad trials
    lo_bin2 = lo_bin2_trials(~ismember(lo_bin2_trials(:, 1), badind), :);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%%% for resampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %response trajectories
    %expected & unexpected (bin1-bin4)
    data.dist_exp1{scnt} = pcat.distance(:, exp_bin1(:, 1));
    data.dist_exp2{scnt} = pcat.distance(:, exp_bin2(:, 1));
    data.dist_un1{scnt} = pcat.distance(:, un_bin1(:, 1));
    data.dist_un2{scnt} = pcat.distance(:, un_bin2(:, 1));
    
    %focused & divided
    data.dist_foc1{scnt} = pcat.distance(:, foc_bin1(:, 1));
    data.dist_foc2{scnt} = pcat.distance(:, foc_bin2(:, 1));
    data.dist_div1{scnt} = pcat.distance(:, div_bin1(:, 1));
    data.dist_div2{scnt} = pcat.distance(:, div_bin2(:, 1));
    
    %high & low coherence
    data.dist_hi1{scnt} = pcat.distance(:, hi_bin1(:, 1));
    data.dist_hi2{scnt} = pcat.distance(:, hi_bin2(:, 1));
    data.dist_lo1{scnt} = pcat.distance(:, lo_bin1(:, 1));
    data.dist_lo2{scnt} = pcat.distance(:, lo_bin2(:, 1));
    
    %response error
    %expected & unexpected (bin1-bin4)
    data.resp_exp1{scnt} = pcat.angerror(:, exp_bin1(:, 1));
    data.resp_exp2{scnt} = pcat.angerror(:, exp_bin2(:, 1));
    data.resp_un1{scnt} = pcat.angerror(:, un_bin1(:, 1));
    data.resp_un2{scnt} = pcat.angerror(:, un_bin2(:, 1));
    
    %focused & divided
    data.resp_foc1{scnt} = pcat.angerror(:, foc_bin1(:, 1));
    data.resp_foc2{scnt} = pcat.angerror(:, foc_bin2(:, 1));
    data.resp_div1{scnt} = pcat.angerror(:, div_bin1(:, 1));
    data.resp_div2{scnt} = pcat.angerror(:, div_bin2(:, 1));
    
    %high & low coherence
    data.resp_hi1{scnt} = pcat.angerror(:, hi_bin1(:, 1));
    data.resp_hi2{scnt} = pcat.angerror(:, hi_bin2(:, 1));
    data.resp_lo1{scnt} = pcat.angerror(:, lo_bin1(:, 1));
    data.resp_lo2{scnt} = pcat.angerror(:, lo_bin2(:, 1));
    
    %interaction effects
    %coh & att
    lo_foc.trials = intersect(lo_bin1(:, 1), foc_bin1(:, 1));
    hi_foc.trials = intersect(hi_bin1(:, 1), foc_bin1(:, 1));
    lo_div.trials = intersect(lo_bin1(:, 1), div_bin1(:, 1));
    hi_div.trials = intersect(hi_bin1(:, 1), div_bin1(:, 1));
    
    %coh & exp
    lo_exp.trials = intersect(lo_bin1(:, 1), exp_bin1(:, 1));
    hi_exp.trials = intersect(hi_bin1(:, 1), exp_bin1(:, 1));
    lo_un.trials = intersect(lo_bin1(:, 1), un_bin1(:, 1));
    hi_un.trials = intersect(hi_bin1(:, 1), un_bin1(:, 1));
    
    %att & exp
    foc_exp.trials = intersect(foc_bin1(:, 1), exp_bin1(:, 1));
    div_exp.trials = intersect(div_bin1(:, 1), exp_bin1(:, 1));
    foc_un.trials = intersect(foc_bin1(:, 1), un_bin1(:, 1));
    div_un.trials = intersect(div_bin1(:, 1), un_bin1(:, 1));
    
    %response trajectories
    %coh & att
    data.dist_lo_foc{scnt} = pcat.distance(:, lo_foc.trials);
    data.dist_hi_foc{scnt} = pcat.distance(:, hi_foc.trials);
    data.dist_lo_div{scnt} = pcat.distance(:, lo_div.trials);
    data.dist_hi_div{scnt} = pcat.distance(:, hi_div.trials);
    %coh & exp
    data.dist_lo_exp{scnt} = pcat.distance(:, lo_exp.trials);
    data.dist_hi_exp{scnt} = pcat.distance(:, hi_exp.trials);
    data.dist_lo_un{scnt} = pcat.distance(:, lo_un.trials);
    data.dist_hi_un{scnt} = pcat.distance(:, hi_un.trials);
    %att & exp
    data.dist_foc_exp{scnt} = pcat.distance(:, foc_exp.trials);
    data.dist_div_exp{scnt} = pcat.distance(:, div_exp.trials);
    data.dist_foc_un{scnt} = pcat.distance(:, foc_un.trials);
    data.dist_div_un{scnt} = pcat.distance(:, div_un.trials);
    
    %response errors
    %coh & att
    data.resp_lo_foc{scnt} = pcat.angerror(:, lo_foc.trials);
    data.resp_hi_foc{scnt} = pcat.angerror(:, hi_foc.trials);
    data.resp_lo_div{scnt} = pcat.angerror(:, lo_div.trials);
    data.resp_hi_div{scnt} = pcat.angerror(:, hi_div.trials);
    %coh & exp
    data.resp_lo_exp{scnt} = pcat.angerror(:, lo_exp.trials);
    data.resp_hi_exp{scnt} = pcat.angerror(:, hi_exp.trials);
    data.resp_lo_un{scnt} = pcat.angerror(:, lo_un.trials);
    data.resp_hi_un{scnt} = pcat.angerror(:, hi_un.trials);
    %att & exp
    data.resp_foc_exp{scnt} = pcat.angerror(:, foc_exp.trials);
    data.resp_div_exp{scnt} = pcat.angerror(:, div_exp.trials);
    data.resp_foc_un{scnt} = pcat.angerror(:, foc_un.trials);
    data.resp_div_un{scnt} = pcat.angerror(:, div_un.trials);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%% plotting stuff by performance bins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %parameters for plotting
    fq = 120; %monitor's refresh rate
    timex = 0: 1000/fq : (frm-1)*(1000/fq); %time axis (x)
    c_foc = [0 0.4 0]; %green
    c_div = [0.4 0 0.6]; %purple
    c_exp = [0 0.4 0.8]; %blue
    c_neu = [0.5 0.5 0.5]; %grey
    c_un = [0.8 0 0]; %orange
    c_hi = [0 0 0]; %black
    c_hi = [0.6 0 0.2]; %maroon
    c_lo = [0.4 0.4 0.4]; %grey
    %         diffcolor1 = [0.5 0.5 0.5]; %exp - unexp
    %         diffcolor2 = [0.5 0.5 0.5];
    %         noise = 0.1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% trajectories plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %as a function of expectation
    for fig = 8*(s-1)+1
        figure(8*(s-1)+2);
        d2 = suptitle(['Response Trajectories: Exp vs. Unexp (subj' num2str(allsub(s)) ')'])
        set(d2, 'Fontsize', 13, 'fontWeight', 'bold')
        set(gca, 'FontSize', 12)
        
        %bin1
        subplot(3, 2, 1) %all exp trials
        plot(timex, pcat.distance(:, exp_bin1(:, 1))); hold on;
        plot(timex, mean(pcat.distance(:, exp_bin1(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
        title('bin1')
        subplot(3, 2, 3) %all unexp trials
        plot(timex, pcat.distance(:, un_bin1(:, 1))); hold on;
        plot(timex, mean(pcat.distance(:, un_bin1(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
        subplot(3, 2, 5) %mean exp vs mean unexp
        plot(timex, mean(pcat.distance(:, exp_bin1(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
        plot(timex, mean(pcat.distance(:, un_bin1(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
        %legend('Expected', 'Unexpected')
        
        %bin2
        if isempty(exp_bin2)
            subplot(3, 2, 2)
            plot(timex, nan);
            title('bin2')
            if isempty(un_bin2)
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, nan);
            else
                subplot(3, 2, 4)
                plot(timex, pcat.distance(:, un_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
            end
        else
            if isempty(un_bin2)
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, exp_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, exp_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, pcat.distance(:, un_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
                plot(timex, mean(pcat.distance(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
                legend('Expected', 'Unexpected')
            end
        end
        
        %using the same ylim for all plots
        for n=1:6
            AX(n) = subplot(3,2,n);
        end
        set(AX, 'YLim', [0 1]); %needs to fix this later, but for now using the max dist of root 2
        set(AX, 'XLim', [0 1500]);
    end
    
    %as a function of attention
    for fig = 8*(s-1)+3
        figure(8*(s-1)+3);
        d3 = suptitle(['Response Trajectories: Foc vs. Div (subj' num2str(allsub(s)) ')'])
        set(d3, 'Fontsize', 13, 'fontWeight', 'bold')
        set(gca, 'FontSize', 12)
        
        %bin1
        subplot(3, 2, 1) %all foc trials
        plot(timex, pcat.distance(:, foc_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.distance(:, foc_bin1(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
        title('bin1')
        subplot(3, 2, 3) %all div trials
        plot(timex, pcat.distance(:, div_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.distance(:, div_bin1(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
        subplot(3, 2, 5) %mean foc vs mean div
        plot(timex, mean(pcat.distance(:, foc_bin1(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
        plot(timex, mean(pcat.distance(:, div_bin1(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
        
        %bin2
        if isempty(foc_bin2)
            subplot(3, 2, 2)
            plot(timex, nan);
            title('bin2')
            if isempty(div_bin2)
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, nan);
            else
                subplot(3, 2, 4)
                plot(timex, pcat.distance(:, div_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
            end
        else
            if isempty(div_bin2)
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, foc_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, foc_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, pcat.distance(:, div_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
                plot(timex, mean(pcat.distance(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
                legend('Focused', 'Divided')
            end
        end
        
        %using the same ylim for all plots
        for n=1:6
            AX(n) = subplot(3,2,n);
        end
        set(AX, 'YLim', [0 1.0]); %needs to fix this later, but for now using the max dist of root 2
        set(AX, 'XLim', [0 1500]);
    end
    
    %as a function of coherence level
    for fig = 8*(s-1)+4
        figure(8*(s-1)+4);
        d4 = suptitle(['Response Trajectories: Hi vs. Lo (subj' num2str(allsub(s)) ')'])
        set(d4, 'Fontsize', 13, 'fontWeight', 'bold')
        set(gca, 'FontSize', 12)
        
        %         for qq = 1:size(lo_bin3, 1)
        %
        %         figure; %hi bin1
        %         %suptitle('hi coh bin1')
        %         plot(timex, pcat.angerror(:, lo_bin3(qq, 1))); hold on;
        %          title(['trial' num2str(qq)]);
        %  end
        %for qq = 1:size(hi_bin1, 1)
        %             if ismember(qq, [1:34])
        %             figure(390);
        %         subplot(2, 17, qq)
        %         plot(timex, pcat.angerror(:, hi_bin1(qq, 1))); hold on;
        %         title(['trial' num2str(qq)]);
        %             elseif ismember(qq, [35:68])
        %                 figure(391);
        %                 subplot(2, 17, qq-34)
        %         plot(timex, pcat.angerror(:, hi_bin1(qq, 1))); hold on;
        %         title(['trial' num2str(qq)]);
        %         elseif ismember(qq, [69:102])
        %                 figure(392);
        %                 subplot(2, 17, qq-68)
        %         plot(timex, pcat.angerror(:, hi_bin1(qq, 1))); hold on;
        %         title(['trial' num2str(qq)]);
        %          elseif ismember(qq, [103:136])
        %                 figure(393);
        %                 subplot(2, 17, qq-102)
        %         plot(timex, pcat.angerror(:, hi_bin1(qq, 1))); hold on;
        %         title(['trial' num2str(qq)]);
        %         elseif ismember(qq, [137:150])
        %                 figure(394);
        %                 subplot(2, 17, qq-136)
        %         plot(timex, pcat.angerror(:, hi_bin1(qq, 1))); hold on;
        %         title(['trial' num2str(qq)]);
        %             end
        %         end
        %%
        %bin1
        subplot(3, 2, 1) %all hi trials
        plot(timex, pcat.distance(:, hi_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.distance(:, hi_bin1(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
        title('bin1')
        subplot(3, 2, 3) %all lo trials
        plot(timex, pcat.distance(:, lo_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.distance(:, lo_bin1(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
        subplot(3, 2, 5) %mean hi vs lo lo
        plot(timex, mean(pcat.distance(:, hi_bin1(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
        plot(timex, mean(pcat.distance(:, lo_bin1(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
        
        %         %bin2
        %         subplot(3, 2, 2)
        %         plot(timex, pcat.distance(:, hi_bin2(1:end, 1))); hold on;
        %         plot(timex, mean(pcat.distance(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
        %         title('bin2')
        %         subplot(3, 2, 4)
        %         plot(timex, pcat.distance(:, lo_bin2(1:end, 1))); hold on;
        %         plot(timex, mean(pcat.distance(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
        %         subplot(3, 2, 6)
        %         plot(timex, mean(pcat.distance(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
        %         plot(timex, mean(pcat.distance(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
        %
        
        %bin2
        if isempty(hi_bin2)
            subplot(3, 2, 2)
            plot(timex, nan);
            title('bin2')
            if isempty(lo_bin2)
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, nan);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, lo_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
                subplot(3, 2, 4)
                plot(timex, mean(pcat.distance(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
            end
        else
            if isempty(lo_bin2)
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, hi_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.distance(:, hi_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, pcat.distance(:, lo_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.distance(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.distance(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
                plot(timex, mean(pcat.distance(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
            end
        end
        
        %using the same ylim for all plots
        for n=1:6
            AX(n) = subplot(3,2,n);
        end
        set(AX, 'YLim', [0 1]); %needs to fix this later, but for now using the max dist of root 2
        set(AX, 'XLim', [0 1500]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% response error plots %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %as a function of expectation
    for fig = 8*(s-1)+5
        figure(8*(s-1)+5);
        %d2 = suptitle(['Response Error: Exp vs. Unexp (subj' num2str(allsub(s)) ')'])
        d2 = suptitle(['Error:Exp/Unexp(subj' num2str(allsub(s)) ')'])
        
        set(d2, 'Fontsize', 13, 'fontWeight', 'bold')
        set(gca, 'FontSize', 12)
        
        %bin1
        subplot(3, 2, 1) %all exp trials
        plot(timex, pcat.angerror(:, exp_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.angerror(:, exp_bin1(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
        %title('bin1')
        subplot(3, 2, 3) %all unexp trials
        plot(timex, pcat.angerror(:, un_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.angerror(:, un_bin1(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
        subplot(3, 2, 5) %mean exp vs mean unexp
        plot(timex, mean(pcat.angerror(:, exp_bin1(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
        plot(timex, mean(pcat.angerror(:, un_bin1(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
        %legend('Expected', 'Unexpected')
        
        %bin2
        if isempty(exp_bin2)
            subplot(3, 2, 2)
            plot(timex, nan);
            %title('bin4')
            if isempty(un_bin2)
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, nan);
            else
                subplot(3, 2, 4)
                plot(timex, pcat.angerror(:, un_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
            end
        else
            if isempty(un_bin2)
                subplot(3, 2, 2)
                plot(timex, pcat.angerror(:, exp_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
                %title('bin4')
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.angerror(:, exp_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
                %title('bin4')
                subplot(3, 2, 4)
                plot(timex, pcat.angerror(:, un_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, exp_bin2(:, 1)), 2), 'color', c_exp, 'LineWidth', 3); hold on;
                plot(timex, mean(pcat.angerror(:, un_bin2(:, 1)), 2), 'color', c_un, 'LineWidth', 3);
            end
        end
        
        %using the same ylim for all plots
        for n=1:6
            AX(n) = subplot(3,2,n);
        end
        set(AX, 'YLim', [0 100]); %needs to fix this later, but for now using the max dist of root 2
        set(AX, 'XLim', [0 1500]);
    end
    
    %as a function of attention
    for fig = 8*(s-1)+6
        figure(8*(s-1)+6);
        d3 = suptitle(['Response Error: Foc vs. Div (subj' num2str(allsub(s)) ')'])
        set(d3, 'Fontsize', 13, 'fontWeight', 'bold')
        set(gca, 'FontSize', 12)
        
        %bin1
        subplot(3, 2, 1) %all foc trials
        plot(timex, pcat.angerror(:, foc_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.angerror(:, foc_bin1(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
        title('bin1')
        subplot(3, 2, 3) %all div trials
        plot(timex, pcat.angerror(:, div_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.angerror(:, div_bin1(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
        subplot(3, 2, 5) %mean foc vs mean div
        plot(timex, mean(pcat.angerror(:, foc_bin1(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
        plot(timex, mean(pcat.angerror(:, div_bin1(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
        
        %bin2
        if isempty(foc_bin2)
            subplot(3, 2, 2)
            plot(timex, nan);
            title('bin2')
            if isempty(div_bin2)
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, nan);
            else
                subplot(3, 2, 4)
                plot(timex, pcat.angerror(:, div_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
            end
        else
            if isempty(div_bin2)
                subplot(3, 2, 2)
                plot(timex, pcat.angerror(:, foc_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.angerror(:, foc_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, pcat.angerror(:, div_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, foc_bin2(:, 1)), 2), 'color', c_foc, 'LineWidth', 3); hold on;
                plot(timex, mean(pcat.angerror(:, div_bin2(:, 1)), 2), 'color', c_div, 'LineWidth', 3);
            end
        end
        
        %using the same ylim for all plots
        for n=1:6
            AX(n) = subplot(3,2,n);
        end
        set(AX, 'YLim', [0 100]); %needs to fix this later, but for now using the max dist of root 2
        set(AX, 'XLim', [0 1500]);
    end
    
    %as a function of coherence level
    for fig = 8*(s-1)+7
        figure(8*(s-1)+7);
        d4 = suptitle(['Response Error: Hi vs. Lo (subj' num2str(allsub(s)) ')'])
        set(d4, 'Fontsize', 13, 'fontWeight', 'bold')
        set(gca, 'FontSize', 12)
        
        %bin1
        subplot(3, 2, 1) %all hi trials
        plot(timex, pcat.angerror(:, hi_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.angerror(:, hi_bin1(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
        title('bin1')
        subplot(3, 2, 3) %all lo trials
        plot(timex, pcat.angerror(:, lo_bin1(1:end, 1))); hold on;
        plot(timex, mean(pcat.angerror(:, lo_bin1(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
        subplot(3, 2, 5) %mean hi vs lo lo
        plot(timex, mean(pcat.angerror(:, hi_bin1(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
        plot(timex, mean(pcat.angerror(:, lo_bin1(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
        
        %bin4
        if isempty(hi_bin2)
            subplot(3, 2, 2)
            plot(timex, nan);
            title('bin2')
            if isempty(lo_bin2)
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, nan);
            else
                subplot(3, 2, 4)
                plot(timex, pcat.angerror(:, lo_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
            end
        else
            if isempty(lo_bin2)
                subplot(3, 2, 2)
                plot(timex, pcat.angerror(:, hi_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, nan);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3);
            else
                subplot(3, 2, 2)
                plot(timex, pcat.angerror(:, hi_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
                title('bin2')
                subplot(3, 2, 4)
                plot(timex, pcat.angerror(:, lo_bin2(:, 1))); hold on;
                plot(timex, mean(pcat.angerror(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
                subplot(3, 2, 6)
                plot(timex, mean(pcat.angerror(:, hi_bin2(:, 1)), 2), 'color', c_hi, 'LineWidth', 3); hold on;
                plot(timex, mean(pcat.angerror(:, lo_bin2(:, 1)), 2), 'color', c_lo, 'LineWidth', 3);
            end
        end
        
        %using the same ylim for all plots
        for n=1:6
            AX(n) = subplot(3,2,n);
        end
        set(AX, 'YLim', [0 120]);
        set(AX, 'XLim', [0 1500]);
    end
    
    %pre-allocate storage
    pcat.distanceRT2locked = nan(51,size(pcat.distance,2));
    pcat.anDiffRT2locked = nan(51,size(pcat.distance,2));
    pcat.anDiffRT1locked = nan(71,size(pcat.distance,2));
    pcat.distanceRT1locked = nan(71,size(pcat.distance,2));
    cd ../..
end

%save('controls_data.mat', 'data', 'frm', 'allsub', 'timex', 'all'); %save out collapsed stuff in pcat
save('patients_data.mat', 'data', 'frm', 'allsub', 'timex', 'all'); %save out collapsed stuff in pcat

