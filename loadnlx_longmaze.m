function [datamtx, csc_mtx, clusters, spike_waveforms, TT_regions] = loadnlx(filename)
%loads data from a neuralynx output folder
%sorted spike data must be renamed with '_sorted' at end of file name

%filename examples:
%D:\2018-04-22_15-48-29
%F:\Adam M\hpc-rsc_lin-track\prelim\m1-6\2018-04-22_15-48-29

%datamtx key
% col#  data
%   1   time
%   2   xpos
%   3   ypos
%   4   head direction (clockwise angle from "north")
%   5   velocity
%   6   acceleration
%   7   stage index (1 bowl, 2 iti, 3 maze, 4 iti, 5 bowl)
%   8   spike event ids (0 is video sample; 1&2 are Left and Right lick)
%   9   RSC ripple power
%   10  HPC ripple power
%   11  Stillness window flag (stillness == 1)
%   12  HPC ripple window number
%   13  Trial number (track stage only)
%   14  RSC ripple window number


%cscmtx key
% col#  data
%   1   time
%   2   rsc csc
%   3   hpc csc
%   4   stage
%   5   stillness flag (1=still)



%% Ensure the existence of necessary files

%load from user-created csv file identifying LFP tt's [RSC_TT,HPC_TT]
LFP_channels = csvread([filename '\csc.csv']);
TT_regions = csvread([filename '\tt_regions.csv']);
place_holds = find(TT_regions==0);
RSC_TTs = TT_regions(place_holds(1):place_holds(2)); RSC_TTs = RSC_TTs(RSC_TTs>0);
HPC_TTs = TT_regions(place_holds(2):place_holds(3)); HPC_TTs = HPC_TTs(HPC_TTs>0);
O_TTs = TT_regions(place_holds(3):place_holds(4)); O_TTs = O_TTs(O_TTs>0);
TT_regions = cell(3,1); TT_regions{1} = RSC_TTs; TT_regions{2} = HPC_TTs; TT_regions{3} = O_TTs;


%% Position and time

% Load position data
[Timestamps, X, Y, Angles] = Nlx2MatVT([filename '\VT1.nvt'], [1 1 1 1 0 0], 0, 1, [] );
datamtx = [Timestamps', X', Y', Angles'];

% Correct time scaling
nlx_start_time = datamtx(1,1);
datamtx(:,1) = fix_time(datamtx(:,1), nlx_start_time);

% Delete false points
datamtx(datamtx(:,2)==0, 2:4) = nan; %dropped points
datamtx(datamtx(:,4)==450, 4) = 90; %neuralynx HD bug
datamtx(ismember(datamtx(:,4), [0 360]), 4) = nan;


%% Flags 

[flag_mtx, stage_idx,  end_sesh_time] = load_flags(filename, datamtx, nlx_start_time);
stage_idx(datamtx(:,1)>=end_sesh_time) = [];
datamtx(datamtx(:,1)>=end_sesh_time,:) = [];    


%% Correct position

%nan behavior outside of stages
iti_idx = ismember(stage_idx, [2 4]);
datamtx(iti_idx, 2:4) = nan;

%correct position during track stage
track_idx = stage_idx == 3;
accepted_area = [0 760 150 300];
datamtx(track_idx, 1:4) = hard_correct_pos(datamtx(track_idx, 1:4),accepted_area(1:2),accepted_area(3:4));

    %if most deleted, set to mean of accepted area (placehold)
    if sum(isnan(datamtx(track_idx, 2)))/length(datamtx(track_idx, 2)) > 0.8
        datamtx(track_idx, 2:3) = ...
            repmat([mean(accepted_area([1 2])) mean(accepted_area([3 4]))],...
            size(datamtx(track_idx, 2:3),1), 1);
    else
        datamtx(track_idx, 1:4) = hard_correct_stempos(datamtx(track_idx, 1:4));
        datamtx(track_idx, 1:4) = resample_jumps(datamtx(track_idx, 1:4), 150, 50);
        datamtx(track_idx, 1:4) = smooth_pos(datamtx(track_idx, 1:4), 20);
        [datamtx(track_idx, 1:4), rotang] = rotate_pos(datamtx(track_idx, 1:4));
    end

%correct position during bowl stages
xdiff_bowl = nan(2,2);
count = 0;

for stage = [1 5]
    count = count+1;
    accepted_area = [250 450 50 150];
    datamtx(stage_idx == stage, 1:4) = hard_correct_pos(datamtx(stage_idx == stage, 1:4),accepted_area(1:2),accepted_area(3:4));
        %if most deleted, set to mean of accepted area (placehold)
        if sum(isnan(datamtx(stage_idx == stage, 2)))/length(datamtx(stage_idx == stage, 2)) > 0.8
            datamtx(stage_idx == stage, 2:3) = ...
                repmat([mean(accepted_area([1 2])) mean(accepted_area([3 4]))],...
                size(datamtx(stage_idx == stage, 2:3),1), 1);
        else
            datamtx(stage_idx == stage, 1:4) = resample_jumps(datamtx(stage_idx == stage, 1:4), 300, 100);
            datamtx(stage_idx == stage, 1:4) = smooth_pos(datamtx(stage_idx == stage, 1:4), 20);
            xdiff_bowl(count, :) = [max(datamtx(stage_idx == stage, 2)) min(datamtx(stage_idx == stage, 2))] ;
            datamtx(stage_idx == stage, 1:4) =  fixed_rotate_pos(datamtx(stage_idx == stage, 1:4), rotang);    
        end
end

%% Trials

%relevant positions
trial_pos = datamtx(stage_idx==3, 1:2);
min_max_pos = [nanmin(trial_pos(:, 2)) nanmax(trial_pos(:, 2))];
length_maze = abs(diff(min_max_pos));

%use position during lick detections to find ends of maze
if sum(flag_mtx(~isnan(flag_mtx(:, 2)), 2)==1)>=2 && sum(flag_mtx(~isnan(flag_mtx(:, 2)), 2)==2)>=2
    
    %compute x position at each lick-detection timestamp
    lick_xpos = interp1(datamtx(:,1), datamtx(:,2), flag_mtx(:,1));
    min_max_x = sort([nanmedian(lick_xpos(flag_mtx(:,2)==1)) nanmedian(lick_xpos(flag_mtx(:,2)==2))]);
    
    %sanity check
    pos_correct = 0;
    if abs(min_max_x(1) - min_max_pos(1)) > length_maze/5
        warning('maze-end position correction')
        pos_correct = 1;
        min_max_x(1) = min_max_pos(1)+length_maze*0.025;
    end
    if abs(min_max_x(2) - min_max_pos(2)) > length_maze/5
        warning('maze-end position correction')
        pos_correct = 1;
        min_max_x(2) = min_max_pos(2)-length_maze*0.025;
    end
    if pos_correct==1
        plot_maze_ends(datamtx(stage_idx==3,:), min_max_x); 
        title(filename)
    end
    
else %or estimate fromm range
    min_max_x = min_max_pos + [length_maze -length_maze].*0.025;
end

x_end_bounds = linspace(min_max_x(1), min_max_x(2), 11);
end_lo = x_end_bounds(2);
end_hi = x_end_bounds(end-1);


%iterate through time points identifying trials
trial_label = nan(size(trial_pos(:,1)));
previous_sect = 0;
cross_times = [];
trial_count = 0;


%find trial boundary crossings
for itime = 1:size(trial_pos,1)
    
    %skip nan
    if isnan(trial_pos(itime,2))
        continue
    end

    %if coming from center
    if ismember(previous_sect, [0 2])       
        current_sect = find_current_sect(trial_pos(itime,2), end_lo, end_hi);
    %if coming from end   
    elseif ismember(previous_sect, [0 1 3])
        current_sect = find_current_sect(trial_pos(itime,2), end_lo, end_hi);
    end
    
    %crossing
    if previous_sect~=0 && previous_sect~=current_sect
        sections = [previous_sect current_sect];
        current_end = sections(ismember(sections, [1 3]));
        cross_times = [cross_times; [trial_pos(itime,1) current_end]];
    end
    
    %shift forward
    previous_sect=current_sect;
end

%number trials 
%based on last crossing on one side before first crossing on the other
while 1
    trial_count = trial_count+1;
    first1 = min(cross_times(cross_times(:,2)==1,1));
    first3 = min(cross_times(cross_times(:,2)==3,1));
    trial_end = max([first1 first3]);
    if first1 < first3
        trial_start = max(cross_times(cross_times(:,2)==1 & cross_times(:,1)<trial_end,1));
        start_side = 1;
    elseif first1 > first3
        trial_start = max(cross_times(cross_times(:,2)==3 & cross_times(:,1)<trial_end,1));
        start_side = 3;
    else
        break
    end
    
    %load
    trial_label(trial_pos(:,1)>=trial_start & trial_pos(:,1)<trial_end) = trial_count;
    
    %delete used crossings
    cross_times = cross_times(cross_times(:,1)>trial_end, :);
    
    %display trials
    %trial_count = trial_count
end

% trial time bounds
all_trial_labels = nan(size(datamtx(:,1)));
all_trial_labels(stage_idx==3) = trial_label;
num_trials = length(unique(all_trial_labels(~isnan(all_trial_labels))));
trial_time_bounds = nan(num_trials,2);
for itrial = 1:num_trials 
    trial_time_bounds(itrial,:) = [nanmin(trial_pos(trial_label==itrial,1)) nanmax(trial_pos(trial_label==itrial,1))];
end
    


%% Normalize position

%find x min and max
trial_minmax_x = nan(num_trials, 2);
for itrial = 1:num_trials
    trial_minmax_x(itrial, :) = [nanmin(datamtx(all_trial_labels==itrial,2))...
        nanmax(datamtx(all_trial_labels==itrial,2))];
end
trial_minmax_x = median(trial_minmax_x);
xdiff_track = trial_minmax_x(2)-trial_minmax_x(1);
xy_ratio = (max(datamtx(stage_idx==3,2))-min(datamtx(stage_idx==3,2)))...
    /(max(datamtx(stage_idx==3,3))-min(datamtx(stage_idx==3,3)));

%normalize x
datamtx(stage_idx==3,2) = datamtx(stage_idx==3,2) - trial_minmax_x(1);
datamtx(stage_idx==3,2) = datamtx(stage_idx==3,2)./xdiff_track;

%proportion y
datamtx(stage_idx==3,3) = (datamtx(stage_idx==3,3) - min(datamtx(stage_idx==3,3)))./...
    (max(datamtx(stage_idx==3,3)) - min(datamtx(stage_idx==3,3)));
datamtx(stage_idx==3,3) = datamtx(stage_idx==3,3)./xy_ratio;
datamtx(stage_idx==3,3) = datamtx(stage_idx==3,3) + ...
    (0.5 - nanmedian(datamtx(stage_idx==3 & datamtx(:,2)>.2 & datamtx(:,2) <.8,3)));

%normalize bowl position
for stage = [1 5]
    
    %ratio
    xy_ratio = (max(datamtx(stage_idx==stage,2))-min(datamtx(stage_idx==stage,2)))...
        /(max(datamtx(stage_idx==stage,3))-min(datamtx(stage_idx==stage,3)));
    if isnan(xy_ratio)
        xy_ratio=1;
    end
    
    %subtract mins from all position points
    datamtx(stage_idx==stage,2) = datamtx(stage_idx==stage,2) - min(datamtx(stage_idx==stage,2));
    datamtx(stage_idx==stage,3) = datamtx(stage_idx==stage,3) - min(datamtx(stage_idx==stage,3));

    %divide all points by their maximum (normalizing to 1)
    max_pos = max(datamtx(stage_idx==stage, [2 3]));
    if any(max_pos==0)
       max_pos = [1 1]; 
    end
    datamtx(stage_idx==stage,2) = datamtx(stage_idx==stage,2)./max_pos(1);
    datamtx(stage_idx==stage,3) = datamtx(stage_idx==stage,3)./max_pos(2);

    %preserve scale
    datamtx(stage_idx==stage,2) = datamtx(stage_idx==stage,2).*xy_ratio;

    %reposition
    xdiff_bowl_hold = max(xdiff_bowl(:,1)) - min(xdiff_bowl(:,2));
    track_bowl_ratio = xdiff_bowl_hold/xdiff_track;
    if track_bowl_ratio==0
        track_bowl_ratio=1;
    end
    datamtx(stage_idx==stage,2) = datamtx(stage_idx==stage,2) - nanmedian(datamtx(stage_idx==stage,2));
    datamtx(stage_idx==stage, 2:3) = datamtx(stage_idx==stage, 2:3).*track_bowl_ratio;
    datamtx(stage_idx==stage, 2:3) = datamtx(stage_idx==stage, 2:3) + [0.5 0.7] - nanmedian(datamtx(stage_idx==stage, 2:3));

end



%% Velocity

velocity_vect = nan(size(datamtx,1),1);
time_jump = 21;
for i = 1:length(datamtx(:,2))-time_jump
    %pos points for evaluating velocity
    p1 = datamtx(i, 2:3);
    p2 = datamtx(i+time_jump, 2:3);

    %time points for evaluating velocity
    t1 = datamtx(i,1);
    t2 = datamtx(i+time_jump,1);
    
    %velocity
    velocity_vect(i+10) = pdist([p1; p2])/(t2-t1);
end
maze_length_in_pixels = 690;
maze_length_in_cm = 91;
pixels_per_cm = maze_length_in_pixels/maze_length_in_cm;
pixels_per_cm_over_time_jump = pixels_per_cm*time_jump;
datamtx(:,5) = velocity_vect;
datamtx(:,5) = datamtx(:,5).*pixels_per_cm_over_time_jump; %cm/s
datamtx(:,5) = datamtx(:,5)./100; %m/s


%% Acceleration

acceleration_vect = nan(size(datamtx,1),1);
for i = 1:length(datamtx(:,2))-(time_jump+1)
    
    %vel points for evaluating velocity
    p1 = velocity_vect(i);
    p2 = velocity_vect(i+time_jump);

    %time points for evaluating velocity
    t1 = datamtx(i,1);
    t2 = datamtx(i+time_jump,1);
    
    %velocity
    acceleration_vect(i+10) = (p2-p1)/(t2-t1);
end
datamtx(:,6) = acceleration_vect;


%% Resample time

datamtx = resample_time(datamtx); %100hz


%% Merge flags

datamtx = sortrows([[datamtx nan(size(datamtx(:,1)))]; [flag_mtx(:,1) nan(size(flag_mtx,1), size(datamtx,2)-1) flag_mtx(:,2)]]);


%% Spike events

%iterate through tt's loading unique clusters
spike_records = [];
spike_waveforms = [];
num_sorted_tts = 0;
for tt_records = 1:8
    
    %load each sorted tt independently
    local_filename = [filename '\TT' num2str(tt_records) '_sorted.ntt'];
    if exist(local_filename, 'file')
        num_sorted_tts = num_sorted_tts+1;
        
        %load nlx
        [TimestampsTT, CellNumbers, Samples] = loadspikes(local_filename, tt_records, nlx_start_time);

        %eliminate spikes occuring after the end of the session
        unq_cells = unique(CellNumbers);
        for ic = 1:length(unq_cells)
            Samples{ic}(:,:,TimestampsTT(CellNumbers==unq_cells(ic))>=end_sesh_time) = [];
        end
        CellNumbers(TimestampsTT>=end_sesh_time) = [];
        TimestampsTT(TimestampsTT>=end_sesh_time) = [];
        
        %soft load times and ids (not sorted)
        spike_records = [spike_records; [TimestampsTT' CellNumbers']];
        spike_waveforms = [spike_waveforms Samples];
    end
end

%combine spike and pos data (if spike data exists)
if num_sorted_tts > 0
    datamtx = sortrows([[datamtx zeros(size(datamtx(:,end)))];...
        [spike_records(:,1) nan(size(spike_records,1), size(datamtx,2)-1) spike_records(:,2)]]);
else
    datamtx = sortrows([datamtx zeros(size(datamtx(:,end)))]);
end
    
% Clusters
if num_sorted_tts > 0
    clusters = unique(datamtx(datamtx(:,end)>0 & ~isnan(datamtx(:,end)),end));
else
    clusters = [];
end


%% Add lick times
flag_rows = nan(size(flag_mtx, 1), size(datamtx, 2));
flag_rows(:,[1 8]) = flag_mtx;
datamtx = sortrows([datamtx; flag_rows]);


%% Interp behavior at spikes and flags
manual_flag_mtx = flag_mtx(flag_mtx(:,2)==0,:);
parse_count = 0;
stage_idx = zeros(size(datamtx(:,1)));
for i = 1:size(manual_flag_mtx,1)+1
    parse_count = parse_count+1;
    
    %after last flag (second bowl)
    if i == size(manual_flag_mtx,1)+1
        idx = datamtx(:,1)>manual_flag_mtx(end,1);
        datamtx(idx,:) = interp_at_nans(datamtx(idx,:));
        stage_idx(idx) = i;
        break
    end
    
    %interp each stage independently
    idx = datamtx(:,1)<=manual_flag_mtx(i,1) & stage_idx==0;
    stage_idx(idx) = i;  
    if sum(~isnan(datamtx(idx,2))) == 0
        continue
    end   
    datamtx(idx,:) = interp_at_nans(datamtx(idx,:));
end

%% Set stage index

%stage index is created above
datamtx(:,7)=stage_idx;


%% Load CSC data

%csc downsample
LFP_downsample = 20; %1 for none

%create csc mtx
csc_mtx = load_lfp(filename, LFP_channels, nlx_start_time, LFP_downsample);

%add datamtx columns for instantaneous power (zscored)
RSC_col = 2;
HPC_col = 3;
for col = [RSC_col HPC_col]
    
    %calculate from RAW csc signal power
    [csc_power, ~, times, windowsize] = bandpower_boxcar(csc_mtx, [csc_mtx(1,1) csc_mtx(end,1)], col, [100 250]);
    nnan_idx = ~isnan(times) & ~isnan(csc_power);
    csc_power = csc_power(nnan_idx); times = times(nnan_idx);
    interped_power = interp1(times, csc_power, datamtx(:,1));
    datamtx = [datamtx interped_power];
end



%% Delete ITI behavior

%datamtx(datamtx(:,1)>flag_mtx(1) & datamtx(:,1)<flag_mtx(2), [2:6]) = nan;
%datamtx(datamtx(:,1)>flag_mtx(3) & datamtx(:,1)<flag_mtx(4), [2:6]) = nan;

datamtx(ismember(datamtx(:,7), [2 4]), 2:6) = nan;
datamtx(ismember(datamtx(:,7), [2 4]), 2:6) = nan;



%% Flag periods of stillness

min_duration = 10; %s
max_speed = 0.01; %m/s

%preallocate
stillness_flags = zeros(size(datamtx(:,1)));
vid_samp_times = datamtx(datamtx(:,8)==0,1);
for ivs = 1:length(vid_samp_times)
    %datamtx rows in current time window
    current_idx = datamtx(:,1)>vid_samp_times(ivs) & datamtx(:,1)<vid_samp_times(ivs)+min_duration;
    %if all velocities during time window are below threshold
    if sum(datamtx(current_idx,5)>max_speed)==0
        %flag time window as 'still'
        stillness_flags(current_idx) = 1;
    end
end
datamtx = [datamtx stillness_flags];



%% Add some columns to csc mtx

csc_mtx = [csc_mtx csc_dmidx(csc_mtx, datamtx, datamtx(:,7))]; %stage
csc_mtx = [csc_mtx csc_dmidx(csc_mtx, datamtx, stillness_flags)]; %stilness



%% Flag HPC ripple windows
%power_threshold = 3; %z
power_threshold = 0.02;

%stillness criteria
still_mean = nanmean(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), 10));
still_std = nanstd(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), 10));

%preallocate
%ripple_flags = datamtx(:,10)>(still_mean + power_threshold*still_std);
%ripple_flags_hold_hpc = zeros(size(datamtx(:,10)));
[~, sort_idx] = sort(datamtx(:,10), 'descend');
ripple_flags = zeros(size(datamtx(:,10)));
ripple_flags(sort_idx(1:round(length(datamtx(:,10))*power_threshold))) = 1;
ripple_flags_hold_hpc = zeros(size(datamtx(:,10)));

%flag ripple windows
for ivs = 1:length(datamtx(:,10))
    if ripple_flags(ivs)==1
        %current timewindow idx
        tw_idx = datamtx(:,1)>datamtx(ivs,1)-windowsize/2 & datamtx(:,1)<datamtx(ivs,1)+windowsize/2;
        ripple_flags_hold_hpc(tw_idx) = 1;
    end
end

%count ripple windows
count = 2; 
for ivs = 1:length(ripple_flags_hold_hpc)
    if ripple_flags_hold_hpc(ivs)==1
        reach = 0;
        while ripple_flags_hold_hpc(ivs+reach)==1
            ripple_flags_hold_hpc(ivs+reach) = count;
            reach = reach+1;
            if ivs+reach <= length(ripple_flags_hold_hpc)
                continue
            else
                break
            end
        end
        count = count+1;
    end
end
ripple_flags_hold_hpc(ripple_flags_hold_hpc>0) = ripple_flags_hold_hpc(ripple_flags_hold_hpc>0)-1;
datamtx = [datamtx ripple_flags_hold_hpc];




%% Add trials column

datamtx = [datamtx nan(size(datamtx,1),1)];
for itrial = 1:size(trial_time_bounds,1)
    datamtx(datamtx(:,1)>=trial_time_bounds(itrial,1) & datamtx(:,1)<trial_time_bounds(itrial,2), end) = itrial;
end



%% Flag RSC ripple windows

%stillness criteria
still_mean = nanmean(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), 9));
still_std = nanstd(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), 9));

%preallocate
%ripple_flags_rsc = datamtx(:,10)>(still_mean + power_threshold*still_std);
%ripple_flags_rsc_hold = zeros(size(datamtx(:,9)));
[~, sort_idx] = sort(datamtx(:,9), 'descend');
ripple_flags_rsc = zeros(size(datamtx(:,9)));
ripple_flags_rsc(sort_idx(1:round(length(datamtx(:,9))*power_threshold))) = 1;
ripple_flags_rsc_hold = zeros(size(datamtx(:,9)));

%flag ripple windows
for ivs = 1:length(datamtx(:,9))
    if ripple_flags_rsc(ivs)==1
        %current timewindow idx
        tw_idx = datamtx(:,1)>datamtx(ivs,1)-windowsize/2 & datamtx(:,1)<datamtx(ivs,1)+windowsize/2;
        ripple_flags_rsc_hold(tw_idx) = 1;
    end
end

%count ripple windows
count = 2; 
for ivs = 1:length(ripple_flags_rsc_hold)
    if ripple_flags_rsc_hold(ivs)==1
        reach = 0; 
        while ripple_flags_rsc_hold(ivs+reach)==1
            ripple_flags_rsc_hold(ivs+reach) = count;
            reach = reach+1;
            if ivs+reach <= length(ripple_flags_rsc_hold)
                continue
            else
                break
            end
        end
        count = count+1;
    end
end
ripple_flags_rsc_hold(ripple_flags_rsc_hold>0) = ripple_flags_rsc_hold(ripple_flags_rsc_hold>0)-1;
datamtx = [datamtx ripple_flags_rsc_hold];



%% Load tt histology labels

%load from user-created file identifying tt's as either RSC_TT (1s), HPC_TT
%(2s), or other (0s)
try
    rsc_hpc_tts = csvread([filename '\tt.csv']);
    RSC_TTs = find(rsc_hpc_tts==1);
    HPC_TTs = find(rsc_hpc_tts==2);
    Other_TTs = find(rsc_hpc_tts==0);
catch
end



%% House cleaning
datamtx = datamtx(datamtx(:,1)<end_sesh_time, :);


end




%% Internal functions

function timeXYD = hard_correct_pos(timeXYD, lohiX, lohiY)
%nan any pos values outside of rectangle defined by lohiX and lohiY
outside_rectX_idx = timeXYD(:,2)<lohiX(1) | timeXYD(:,2)>lohiX(2); 
outside_rectY_idx = timeXYD(:,3)<lohiY(1) | timeXYD(:,3)>lohiY(2);

%nan points
timeXYD(outside_rectX_idx | outside_rectY_idx, [2 3]) = nan;

end
function timeXYD = hard_correct_stempos(timeXYD)
%nan positions > 3 STDs from the median y pos of the center 3/5ths of the 
%maze stem

%stem pos
min_max_x = [min(timeXYD(:,2)) max(timeXYD(:,2))];
stem_sects = (abs(diff(min_max_x))/5);
stem_pos_idx = timeXYD(:,2)>(min_max_x(1)+stem_sects) & timeXYD(:,2)<(min_max_x(2)-stem_sects);
stem_pos = timeXYD(stem_pos_idx, :);

%y_limits
med_y = median(stem_pos(:,3));
std_y = std(stem_pos(:,3));
lohiY = [med_y-std_y*2 med_y+std_y*2];

%nan any pos values outside of rectangle defined by lohiX and lohiY
outside_rectY_idx = timeXYD(:,3)<lohiY(1) | timeXYD(:,3)>lohiY(2);

%nan points
timeXYD(outside_rectY_idx & stem_pos_idx, [2 3]) = nan;

end
function [timeXYD, dists] = resample_jumps(timeXYD, too_big_origin, cum_mod)
%delete eroneous spatial positions

    %prep to remove improbable changes in position
    too_big = too_big_origin;
    guilty_by_association = 4;

    %adjacent pos points for evaluating velocity
    first_nnan_row = find(~isnan(timeXYD(:,2)),1,'first');
    p1 = timeXYD(first_nnan_row, 2:3);
    p2 = timeXYD(first_nnan_row+1, 2:3);

    %adjacent time points for evaluating velocity
    t1 = timeXYD(first_nnan_row,1);
    t2 = timeXYD(first_nnan_row+1,1);

    %preallocate
    deletions = zeros(size(timeXYD,1),1);
    dists = zeros(size(timeXYD,1),1);

    %iterate through adjacent points to evaluate velocity
    count = 0;
    for i = first_nnan_row:length(timeXYD(:,2))-2

        %velocity
        current_distance = pdist([p1; p2])/(t2-t1);
        dists(i+1) = current_distance;

        %if the current velocity is too big
        if ~(current_distance <= too_big)

            %note that point (and the next 4) should be deleted (index for later)
            if length(timeXYD(:,2))-2-i > guilty_by_association
                deletions(i:i+guilty_by_association) = 1;
            end

            %move to the next point, but keep the first of the adjacent pair
            p2 = timeXYD(i+2, 2:3);
            t2 = timeXYD(i+2, 1);

            %each time it's too big, increase what is considered "too big"
            count = count + cum_mod;
            too_big = too_big + count;

        %if it's not too big
        else

            %reset what is considered "too big"
            too_big = too_big_origin;
            count = 0;

            %update points
            p1 = timeXYD(i+1, 2:3);
            p2 = timeXYD(i+2, 2:3);
            t1 = timeXYD(i+1, 1);
            t2 = timeXYD(i+2, 1);

        end
    end

    
    %figure; histogram(dists,100)
    
    %index to delete dubious points
    deletions = logical(deletions);
    
    %{
    timeidx = timeXYD(:,1)>1200 & timeXYD(:,1)<1400;
    figure; plot3(timeXYD(timeidx,2), timeXYD(timeidx,3), timeXYD(timeidx,1)); title orig
        hold on; plot3(timeXYD(timeidx & ~deletions,2), timeXYD(timeidx & ~deletions,3), timeXYD(timeidx & ~deletions,1), 'b.');
    hold on; plot3(timeXYD(timeidx & deletions,2), timeXYD(timeidx & deletions,3), timeXYD(timeidx & deletions,1), 'r.');
    %}
    
    timeXYD(deletions, 2:3) = NaN;

    %replace deleted points with interpolated values
    non_nan_pos = ~isnan(timeXYD(:,2)) & ~isnan(timeXYD(:,3)); %index
    timeXYD(:,2) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 2), timeXYD(:,1), 'linear');
    timeXYD(:,3) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 3), timeXYD(:,1), 'linear');
    
    %figure; plot3(timeXYD(timeidx,2), timeXYD(timeidx,3), timeXYD(timeidx,1)); title after

end
function timeXYD = resample_time(timeXYD)
    %resample time and interpolate new position values
    pos(:,1) = (0:0.01:max(timeXYD(:,1)))';
    pos(:,2) = interp1(timeXYD(:,1), timeXYD(:,2), pos(:,1), 'linear');
    pos(:,3) = interp1(timeXYD(:,1), timeXYD(:,3), pos(:,1), 'linear');
    pos(:,4) = interp_hd(timeXYD(:,1), timeXYD(:,4), pos(:,1), 'linear');
    pos(:,5) = interp1(timeXYD(:,1), timeXYD(:,5), pos(:,1), 'linear');
    pos(:,6) = interp1(timeXYD(:,1), timeXYD(:,6), pos(:,1), 'linear');
    timeXYD = pos;
end
function [TimestampsTT, CellNumbers, Waveforms] = loadspikes(filename, tt_records, nlx_start_time)
%use nlx functions to load spike data from a single sorted cluster file

    %nlx load function
    [TimestampsTT, CellNumbers, Samples] = Nlx2MatSpike(filename, ...
        [1 0 1 0 1], 0, 1, [] );

    %delete noise
    Samples(:,:,CellNumbers==0) = [];
    TimestampsTT(CellNumbers==0) = [];
    CellNumbers(CellNumbers==0) = [];
    
    %clusters
    cluster_numbers = unique(CellNumbers);

    %index and load samples
    %change id for each neuron to tt.cluster
    Waveforms = cell(1, length(cluster_numbers));
    unique_cell = 0;
    for ic = cluster_numbers
        unique_cell = unique_cell+1;
        Waveforms{unique_cell} = Samples(:,:,CellNumbers==ic);
        CellNumbers(CellNumbers==ic) = tt_records + unique_cell/100;
    end
    
    %correct time
    TimestampsTT = fix_time(TimestampsTT, nlx_start_time);

end
function fixed_timestamps = fix_time(timestamps, start_time)
    timestamps = timestamps - start_time;
    fixed_timestamps = timestamps./1000000;
end
function interped_hd = interp_hd(times_orig, HDs, times_new, varargin)
% INTERP_LON interpolates a set of longitude angles (in deg)
%
% Usage: out = interp_lon(x,lon,xq)
%
% x and lon are vectors of length N.  function evalutes longitude 
% (in deg -180..180) at points xq using unwrap and interp1
%
% to specify interpolation method used in interp1, use
% out = interp_lon(x,lon,xq,METHOD)
%

    HDs = HDs - 180; %ampm
    
    ulon=unwrap(HDs*pi/180)*180/pi;
    if nargin>3
      interped_hd=interp1(times_orig,ulon,times_new,varargin{1});
    else
      interped_hd=interp1(times_orig,ulon,times_new);
    end
    interped_hd=mod(interped_hd,360);
    interped_hd(interped_hd>180)=interped_hd(interped_hd>180)-360;
    
    interped_hd = interped_hd + 180; %ampm
end
function smoothed_hd = smooth_hd(HDs, smoothing_window_size)
%based on interp_hd
    HDs = HDs - 180;
    ulon=unwrap(HDs*pi/180)*180/pi;
    smoothed_hd=smooth(ulon,smoothing_window_size);
    smoothed_hd=mod(smoothed_hd,360);
    smoothed_hd(smoothed_hd>180)=smoothed_hd(smoothed_hd>180)-360;
    smoothed_hd = smoothed_hd + 180;
end
function timeXYD = interp_at_nans(timeXYD)
%replace deleted points with interpolated values
    non_nan_pos = ~isnan(timeXYD(:,2)) & ~isnan(timeXYD(:,3)); %pos index
    non_nan_hd = ~isnan(timeXYD(:,4)); %hd index
    timeXYD(:,2) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 2), timeXYD(:,1), 'linear');
    timeXYD(:,3) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 3), timeXYD(:,1), 'linear');
    
    if ~isempty(timeXYD(non_nan_hd, 1))
        timeXYD(:,4) = interp_hd(timeXYD(non_nan_hd, 1), timeXYD(non_nan_hd, 4), timeXYD(:,1), 'linear');
    end
    
    timeXYD(:,5) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 5), timeXYD(:,1), 'linear');
    timeXYD(:,6) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 6), timeXYD(:,1), 'linear');
end
function timeXYD = smooth_pos(timeXYD, smoothing_window_size)
    non_nan_pos = ~isnan(timeXYD(:,2)) & ~isnan(timeXYD(:,3)) & ~isnan(timeXYD(:,4)); %index
    timeXYD(non_nan_pos,2) = smooth(timeXYD(non_nan_pos, 2), smoothing_window_size);
    timeXYD(non_nan_pos,3) = smooth(timeXYD(non_nan_pos, 3), smoothing_window_size);
    timeXYD(non_nan_pos,4) = smooth_hd(timeXYD(non_nan_pos, 4), smoothing_window_size');
end
function [datamtx, rotang] = rotate_pos(datamtx)
%rotates positions and and HDs so the track (mouse's path) is parallel with
%the x axis

%axis bins
bins = 50;
xaxis_bins = linspace(floor(min(datamtx(:,2))),ceil(max(datamtx(:,2))),bins);
xy_ratio = (max(datamtx(:,2))-min(datamtx(:,2)))/(max(datamtx(:,3))-min(datamtx(:,3)));

%median positions
median_x = nan(2,bins-1);

    %for each bin (in both axes)
    for bin = 1:(bins-1)

        %axis indices
        x_index = datamtx(:,2)>xaxis_bins(bin) & datamtx(:,2)<xaxis_bins(bin+1);

        %find median for x and y
        median_x(2,bin) = nanmedian(datamtx(x_index, 3));
    end 
    median_x(1,:) = xaxis_bins(1:end-1) + (xaxis_bins(2)-xaxis_bins(1))/2;

    %negative slope of line fit to position points
    angle_offset = polyfit(median_x(1,:), median_x(2,:), 1);
    angle_offset = angle_offset(1);

    %ROTATION
    %corrects for the misalignment of the track
    %

    %counterclockwise rotation angle in radians
    rotang = -(atan(angle_offset));

    %build rotation matrix
    romat=[cos(rotang) -sin(rotang);sin(rotang) cos(rotang)];       

    %rotation occurs around origin, so we temporarily center all points around the origin
    datamtx(:,2) = datamtx(:,2) - mean(median_x(1,:));
    datamtx(:,3) = datamtx(:,3) - mean(median_x(2,:));

    %apply rotation
    datamtx(:,[2 3]) = (romat*datamtx(:,[2 3])')';

    %undo centering
    datamtx(:,2) = datamtx(:,2) + mean(median_x(1,:));
    datamtx(:,3) = datamtx(:,3) + mean(median_x(2,:));      
    
    %interp deleted points
    %  
    %index for missing position values
    real_x = datamtx(~isnan(datamtx(:,2)),2);
    real_y = datamtx(~isnan(datamtx(:,3)),3);
    all_time = datamtx(:,1);

    %only including unique time points among the missing rows
    [real_time_x, idx_x] = unique(datamtx(~isnan(datamtx(:,2)),1));
    [real_time_y, idx_y] = unique(datamtx(~isnan(datamtx(:,3)),1));

    %index real_x and real_y to match real_time unique elements
    datamtx(:,2) = interp1(real_time_x, real_x(idx_x), all_time);
    datamtx(:,3) = interp1(real_time_y, real_y(idx_y), all_time); 
    
    %rotate HDs
    datamtx(:,4) = datamtx(:,4) + rad2deg(rotang);
    
end
function [datamtx] = fixed_rotate_pos(datamtx, rotang)
    %xy dimension ranges and ratio
    xy_ratio = (max(datamtx(:,2))-min(datamtx(:,2)))/(max(datamtx(:,3))-min(datamtx(:,3)));

    %build rotation matrix
    romat=[cos(rotang) -sin(rotang);sin(rotang) cos(rotang)];       

    %rotation occurs around origin, so we temporarily center all points around the origin
    median_xy = nanmedian(datamtx(:,2:3));
    datamtx(:,2) = datamtx(:,2) - median_xy(1);
    datamtx(:,3) = datamtx(:,3) - median_xy(2);
    
    %apply rotation
    datamtx(:,[2 3]) = (romat*datamtx(:,[2 3])')';

    %undo centering
    datamtx(:,2) = datamtx(:,2) + median_xy(1);
    datamtx(:,3) = datamtx(:,3) + median_xy(2);      
            
    %SHIFT x positions to [new_xdims] and y positions to vary around ycenter
    %
    center_y = median_xy(2);

    %interp deleted points
    %  
    %index for missing position values
    real_x = datamtx(~isnan(datamtx(:,2)),2);
    real_y = datamtx(~isnan(datamtx(:,3)),3);
    all_time = datamtx(:,1);

    %only including unique time points among the missing rows
    [real_time_x, idx_x] = unique(datamtx(~isnan(datamtx(:,2)),1));
    [real_time_y, idx_y] = unique(datamtx(~isnan(datamtx(:,3)),1));
    
    %index real_x and real_y to match real_time unique elements
    datamtx(:,2) = interp1(real_time_x, real_x(idx_x), all_time);
    datamtx(:,3) = interp1(real_time_y, real_y(idx_y), all_time); 
    
    %rotate HDs
    datamtx(:,4) = datamtx(:,4) + rad2deg(rotang);
    
end
function csc_mtx = load_lfp(filename, LFP_channels, nlx_start_time, downsample)
%loads csc sample data from CSC channels listed in LFP_channels
%csc_mtx is a column of timestamps followed by 1 column for each channel
%keep one CSC sample out of every number indicated by 'downsample'

%preallocate
csc_mtx = [];

    %iterate through channels
    count = 0;
    for ilfp = LFP_channels
        count = count+1;

        %filename
        fn = [filename '\CSC' num2str(ilfp) '.ncs']

        %neuralynx reader
        [Timestamps, Samples] = Nlx2MatCSC(fn,[1 0 0 0 1], 0, 1, [] );

        %on first iteration
        if count == 1
            %correct timestamps
            Timestamps = (Timestamps-nlx_start_time)./1000000;
            
            %interpolate timestamps to full sample rate
            Timestamps = linspace(min(Timestamps), max(Timestamps), length(Timestamps)*512);
            
            %downsample
            Timestamps = Timestamps(1:downsample:end);
            
            %load timestamps
            csc_mtx = [csc_mtx Timestamps'];
            
            %preallocation
            csc_mtx = [csc_mtx zeros(size(csc_mtx,1), length(LFP_channels))];
        end

        %load downsampled samples
        Samples = Samples(:);
        csc_mtx(:,count+1) = Samples(1:downsample:end);
    end
end
function [flag_mtx, stage_idx, end_sesh_time] = load_flags(filename, datamtx, nlx_start_time)
% integrate flag input into datamtx

    %default end sesh at last time points
    end_sesh_time = datamtx(end,1);

    %matlab reader
    [Timestamps_events, EventIDs, EventStrings] = Nlx2MatEV([filename '\Events.nev'], [1 1 0 0 1], 0, 1, [] );
    
    %delete neuralynx automatic events
    Timestamps_events(EventIDs==19) = []; 
    EventStrings(EventIDs==19) = [];
    EventIDs(EventIDs==19) = [];

    %EventString index
    EvStrIdx = nan(size(EventIDs));
    EvStrIdx(EventIDs==4) = 0; %manual flag
    EvStrIdx(contains(EventStrings, '01).')) = 1; %first ttl (left rwd)
    EvStrIdx(contains(EventStrings, '02).')) = 2; %second ttl (right rwd)
    
    %delete extra lick detection flags
    Timestamps_events(isnan(EvStrIdx)) = []; 
    EventIDs(isnan(EvStrIdx)) = [];
    EvStrIdx(isnan(EvStrIdx)) = [];
    
    %correct time
    Timestamps_events = (Timestamps_events - nlx_start_time)./1000000;
    
    %merge time and flag index
    flag_mtx = [Timestamps_events' EvStrIdx'];
    
    %check for correct tnumber of manual flags
    if sum(EventIDs==4)==8
        display([num2str(sum(EventIDs==4)) ' manual events found. Excluding data after 5th flag.'])
        end_sesh_time = max(Timestamps_events(find(EvStrIdx==0,5)))
        flag_mtx(flag_mtx(:,1)>=max(Timestamps_events(find(EvStrIdx==0,5))),:) = []
    elseif sum(EventIDs==4)~=4
        display([num2str(sum(EventIDs==4)) ' manual events found. Only using the greatest spaced 4.'])
    end

    %parse datamtx using 4 manual flags (end bowl, start maze, end maze,
    %start bowl)
    stage_idx = nan(size(datamtx,1), 1);
    manual_flag_times = Timestamps_events(EventIDs==4);
    %manual_flag_times(2) = [];
    
    %TEST ampm 1/9/2019
    if length(manual_flag_times)~=4
        figure; plot(manual_flag_times, 1, 'o')
        length(manual_flag_times)
       spacing = manual_flag_times(2:end) - manual_flag_times(1:end-1);
       manual_flag_times(find(spacing==max(spacing))+1) = [];
    end
    
    for ievent = 1:4
       stage_idx(datamtx(:,1)<= manual_flag_times(ievent) & isnan(stage_idx)) = ievent;
       if ievent == 4
           stage_idx(datamtx(:,1) > manual_flag_times(ievent)) = ievent+1;
       end
    end
end
function plot_maze_ends(datamtx, maze_ends)
    figure; hold on    
    plot(datamtx(:,2), datamtx(:,3))
    plot(maze_ends(1), mean(datamtx(:,3)), 'r.', 'markersize', 25)
    plot(maze_ends(2), mean(datamtx(:,3)), 'r.', 'markersize', 25)
end
function current_sect = find_current_sect(trial_pos_x, end_lo, end_hi)
    %if on the center
    if trial_pos_x>end_lo && trial_pos_x<end_hi
        current_sect = 2;
    %if on L
    elseif trial_pos_x<end_lo
        current_sect = 1;
    %if on R
    elseif trial_pos_x>end_hi
        current_sect = 3;
    end
end
