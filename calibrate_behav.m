% This script calibrates subjects' joystick responses based on their
% averaged representation of each stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
scnt = 0;
numsess = 2; %2 sessions/subj
frm = 500;
%sub = [11 12 13 14 15 16 17 18 19 20 21 24]; % controls
sub = [30 31 32 33] % patients

for s = 1:size(sub, 2) %subject #
    s
    scnt = scnt+1;
    clear pcat in_file out_file
    
    cd (['data/subj' num2str(sub(s))]);
    
    if sub(s) < 10
        in_file = dir(['AxeHC_Subj0' num2str(sub(s)) '*Cali1*']);
        %set up a name for an output file
        out_file = ['AxeHC_Cali_Subj0' num2str(sub(s)) '.mat'];
    elseif sub(s) > 9
        in_file = dir(['AxeHC_Subj' num2str(sub(s)) '*Cali1*']);
        
        %output file
        out_file = ['AxeHC_Cali_Subj' num2str(sub(s)) '.mat'];
    end
    
    label = {'stimDir', 'stimDirReal', 'joyx', 'joyy'};
    %stimDir = tag assigned to each direction (1-7)
    %stimDirReal = actual direction in rads
    %joyx/joyy = x/y coordinate of the joystick
    
    %loop through all runs
    for r = 1:size(in_file,1)
        clear p
        load(in_file(r).name);
        
        %loop through all 60 trials of calibration
        for t = 1:60
            %distance of the joystick movement
            p.joyxx(:,t) = p.joyx(1:1000, t);
            p.joyyy(:,t) = p.joyy(1:1000, t);
            
            dd(:, t) = sqrt(p.joyxx(:,t).^2 + p.joyyy(:,t).^2);
            %find index for max distance
            [~,ind(t)] = max(dd(:, t));
            
            %find (x,y) at max distance
            x_raw(t) = p.joyxx(ind(t), t);
            y_raw(t) = p.joyyy(ind(t), t);
            
            %
            distc = 1; %corrected joystick distance: 1 au
            rampup_ind = 1:ind(t);
            [c index(t)] = min(abs(dd(rampup_ind, t) - distc)); %use this instead to make sure it happens before the peak
            
            %find (x,y) at max corrected distance which is ~ 1 au
            x(1, t) = p.joyxx(index(t), t);
            y(1, t) = p.joyyy(index(t), t);
            
            if index(t) <= frm %response was made before deadline
                %111918
                %if ind(t) <= frm; %try defining late responses in a different way
                x(2, t) = p.joyxx(index(t), t); %2nd row says nan for late responses
                y(2, t) = p.joyyy(index(t), t);
                
            else %late response
                x(2, t) = NaN;
                y(2, t) = NaN;
            end
            %
            %correcting for the joysitck's square base
            %             distc = 1; %corrected joystick distance: 1 au
            %             [c index] = min(abs(d - distc));
            %             %find (x,y) at max corrected distance which is ~ 1 au
            %             x(t) = p.joyx(index, t);
            %             y(t) = p.joyy(index, t);
            
        end
        %find angle of joystick trajectory based on (x,y) at max distance
        [angles, disp] = (cart2pol(x,y));
        ang(r, :) = wrapToPi(angles(1, :)); %need to wrap to [-pi pi] so it can later be fed to circ_mean
        
        %concatenate
        if r == 1
            for ii = 1:2
                pcat.(label{ii}) = p.(label{ii});
            end
            %concatenate angles of joystick trajectory
            pcat.angle = ang(r, :);
            
        else
            for ii = 1:2
                pcat.(label{ii}) = cat(2, pcat.(label{ii}), p.(label{ii}));
            end
            pcat.angle = cat(2, pcat.angle, ang(r, :));
        end
    end
    
    %loop through all 5 directions
    for j = 1:5
        cali.meanangle(j) = circ_mean(pcat.angle(pcat.stimDir==j)'); %mean
        cali.medianangle(j) = circ_median(pcat.angle(pcat.stimDir==j)'); %median
        cali.stdangle(j) = circ_std(pcat.angle(pcat.stimDir==j)'); %std
    end
    
    %% plot parameters
    cmap = colormap(parula(5));
    
    %% plot all trials
    %close all
    figure(s);
    %suptitle('poop')
    suptitle(['\fontsize{20} Calibration (Subj', num2str(sub(s)) ')'])
    
    subplot(1,3,1)
    polar(0, 4, 'w')
    hold on
    
    %response directions (for all sessions)
    A1 = polar(pcat.angle(pcat.stimDir==1), repmat(0.5, 1, numsess*12) , 'o'); hold on;
    A2 = polar(pcat.angle(pcat.stimDir==2), repmat(1, 1, numsess*12), 'o'); hold on;
    A3 = polar(pcat.angle(pcat.stimDir==3), repmat(1.5, 1, numsess*12), 'o'); hold on;
    A4 = polar(pcat.angle(pcat.stimDir==4), repmat(2, 1, numsess*12), 'o'); hold on;
    A5 = polar(pcat.angle(pcat.stimDir==5), repmat(2.5, 1, numsess*12), 'o'); hold on;
    
    set(A1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
    set(A2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
    set(A3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
    set(A4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
    set(A5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))
    
    %actual presented directions
    B1 = polar(p.stimDirReal(1), 0.5, 'd'); hold on;
    B2 = polar(p.stimDirReal(2), 1, 'd'); hold on;
    B3 = polar(p.stimDirReal(3), 1.5, 'd'); hold on;
    B4 = polar(p.stimDirReal(4), 2, 'd'); hold on;
    B5 = polar(p.stimDirReal(5), 2.5, 'd'); hold on;
    
    set(B1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
    set(B2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
    set(B3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
    set(B4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
    set(B5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))
    
    hold off
    title(['\fontsize{16} All trials'])
    
    
    %% plot median
    subplot(1,3,2)
    polar(0, 4, '-k')
    hold on
    
    %response directions
    A1 = polar(cali.medianangle(1), 0.5, 'o'); hold on;
    A2 = polar(cali.medianangle(2), 1, 'o'); hold on;
    A3 = polar(cali.medianangle(3), 1.5, 'o'); hold on;
    A4 = polar(cali.medianangle(4), 2, 'o'); hold on;
    A5 = polar(cali.medianangle(5), 2.5, 'o'); hold on;
    
    set(A1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
    set(A2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
    set(A3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
    set(A4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
    set(A5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))
    
    %actual presented directions
    B1 = polar(p.stimDirReal(1), 0.5, 'd'); hold on;
    B2 = polar(p.stimDirReal(2), 1, 'd'); hold on;
    B3 = polar(p.stimDirReal(3), 1.5, 'd'); hold on;
    B4 = polar(p.stimDirReal(4), 2, 'd'); hold on;
    B5 = polar(p.stimDirReal(5), 2.5, 'd'); hold on;
    
    set(B1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
    set(B2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
    set(B3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
    set(B4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
    set(B5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))
    
    title(['\fontsize{16} Median'])
    
    %% plot mean
    subplot(1,3,3)
    polar(0, 4, '-k')
    hold on
    
    %response directions
    A1 = polar(cali.meanangle(1), 0.5, 'o'); hold on;
    A2 = polar(cali.meanangle(2), 1, 'o'); hold on;
    A3 = polar(cali.meanangle(3), 1.5, 'o'); hold on;
    A4 = polar(cali.meanangle(4), 2, 'o'); hold on;
    A5 = polar(cali.meanangle(5), 2.5, 'o'); hold on;
    
    set(A1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
    set(A2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
    set(A3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
    set(A4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
    set(A5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))
    
    %actual presented directions
    B1 = polar(p.stimDirReal(1), 0.5, 'd'); hold on;
    B2 = polar(p.stimDirReal(2), 1, 'd'); hold on;
    B3 = polar(p.stimDirReal(3), 1.5, 'd'); hold on;
    B4 = polar(p.stimDirReal(4), 2, 'd'); hold on;
    B5 = polar(p.stimDirReal(5), 2.5, 'd'); hold on;
    
    set(B1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
    set(B2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
    set(B3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
    set(B4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
    set(B5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))
    
    title(['\fontsize{16} Mean'])
    
    %%
    for iii = 1:5
        allmean(scnt, iii) = cali.meanangle(:, iii);
        allmedian(scnt, iii) = cali.medianangle(:, iii);
        allstd(scnt, iii) = cali.stdangle(:, iii);
    end
    
    
    %% save the output file
    save([out_file '.mat'], 'cali')
    cd ../..
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% plot across-subject averages %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allmean2 = wrapToPi(allmean);
allmedian2 = wrapToPi(allmedian);
allstd2 = wrapToPi(allstd);

avgmean = circ_mean(allmean2, [], 1);
avgmedian = circ_median(allmedian2);
%% plot median
figure(size(sub, 2)+1)
subplot(1,2,1)
polar(0, 4, '-k')
hold on

%response directions
A1 = polar(avgmedian(1), 0.5, 'o'); hold on;
A2 = polar(avgmedian(2), 1, 'o'); hold on;
A3 = polar(avgmedian(3), 1.5, 'o'); hold on;
A4 = polar(avgmedian(4), 2, 'o'); hold on;
A5 = polar(avgmedian(5), 2.5, 'o'); hold on;

set(A1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
set(A2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
set(A3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
set(A4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
set(A5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))

%actual presented directions
B1 = polar(p.stimDirReal(1), 0.5, 'd'); hold on;
B2 = polar(p.stimDirReal(2), 1, 'd'); hold on;
B3 = polar(p.stimDirReal(3), 1.5, 'd'); hold on;
B4 = polar(p.stimDirReal(4), 2, 'd'); hold on;
B5 = polar(p.stimDirReal(5), 2.5, 'd'); hold on;

set(B1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
set(B2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
set(B3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
set(B4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
set(B5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))

title(['\fontsize{16} Median'])

%% plot mean
subplot(1,2,2)
polar(0, 4, '-k')
hold on

%response directions
A1 = polar(avgmean(1), 0.5, 'o'); hold on;
A2 = polar(avgmean(2), 1, 'o'); hold on;
A3 = polar(avgmean(3), 1.5, 'o'); hold on;
A4 = polar(avgmean(4), 2, 'o'); hold on;
A5 = polar(avgmean(5), 2.5, 'o'); hold on;

set(A1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
set(A2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
set(A3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
set(A4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
set(A5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))

%actual presented directions
B1 = polar(p.stimDirReal(1), 0.5, 'd'); hold on;
B2 = polar(p.stimDirReal(2), 1, 'd'); hold on;
B3 = polar(p.stimDirReal(3), 1.5, 'd'); hold on;
B4 = polar(p.stimDirReal(4), 2, 'd'); hold on;
B5 = polar(p.stimDirReal(5), 2.5, 'd'); hold on;

set(B1, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(1, :))
set(B2, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(2, :))
set(B3, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(3, :))
set(B4, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(4, :))
set(B5, 'MarkerSize', 8, 'LineWidth',3, 'color', cmap(5, :))

title(['\fontsize{16} Mean'])

