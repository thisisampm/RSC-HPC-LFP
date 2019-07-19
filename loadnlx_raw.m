function [datamtx_txyhs, frm_rng, datamtx_all] = loadnlx_raw(folderpath, varargin)
%% COMPUTE raw times, xy_position, head direction, and spike times

if nargin > 1
    time_range = varargin{1};
end

%% Position

% Load position data
[Timestamps, X, Y, Angles] = Nlx2MatVT([folderpath '\VT1.nvt'], [1 1 1 1 0 0], 0, 1, [] );
datamtx = [Timestamps', X', Y', Angles'];

% Correct time scaling
nlx_start_time = datamtx(1,1);
datamtx(:,1) = fix_time(datamtx(:,1), nlx_start_time);

%time range index
time_range_idx = datamtx(:,1)>=time_range(1) & datamtx(:,1)<=time_range(2);


% Delete false points
datamtx(datamtx(:,2)==0, 2:4) = nan; %dropped points
datamtx(datamtx(:,4)==450, 4) = 90; %neuralynx HD bug
datamtx(ismember(datamtx(:,4), [0 360]), 4) = nan;

%flags
[flag_mtx, stage_idx,  end_sesh_time] = load_flags(folderpath, datamtx, nlx_start_time);
stage_idx(datamtx(:,1)>=end_sesh_time) = nan;
datamtx(datamtx(:,1)>=end_sesh_time,:) = nan;    


%correct position
%
%nan behavior outside of stages
iti_idx = ismember(stage_idx, [2 4]);
datamtx(iti_idx, 2:4) = nan;


if any(stage_idx(time_range_idx)==3)

    %correct position during track stage
    track_idx = stage_idx == 3;
    accepted_area = [75 650 275 365];
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
        end
    
end

%correct position during bowl stages
xdiff_bowl = nan(2,2);
count = 0;

for stage = [1 5]

    if ~any(stage_idx(time_range_idx)==stage)
        continue
    end
    
    count = count+1;
    accepted_area = [125 450 125 250];
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
            %datamtx(stage_idx == stage, 1:4) =  fixed_rotate_pos(datamtx(stage_idx == stage, 1:4), rotang);    
        end
end


%% COMPUTE Spike event times

%iterate through tt's loading unique clusters
spike_records = [];
spike_waveforms = [];
num_sorted_tts = 0;
for tt_records = 1:8
    
    %load each sorted tt independently
    local_filename = [folderpath '\TT' num2str(tt_records) '_sorted.ntt'];
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


%% Trials

%relevant positions
trial_pos = datamtx(stage_idx==3, 1:2);
min_max_pos = [nanmin(trial_pos(:, 2)) nanmax(trial_pos(:, 2))];
length_maze = abs(diff(min_max_pos));

%use position during lick detections to find ends of maze
if sum(flag_mtx(~isnan(flag_mtx(:, 2)), 2)==1)>=2 && sum(flag_mtx(~isnan(flag_mtx(:, 2)), 2)==2)>=2    
    
    
    %compute x position at each lick-detection timestamp
    [unq1, unq1_idx] = unique(datamtx(:,1));
    unq1_idx = ismember(1:size(datamtx,1), unq1_idx)';
    nnan_idx = ~isnan(datamtx(:,1)) & ~isnan(datamtx(:,2));
    
    lick_xpos = interp1(datamtx(nnan_idx & unq1_idx,1), datamtx(nnan_idx & unq1_idx,2), flag_mtx(:,1));
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
    trial_count = trial_count;
end

% trial time bounds
all_trial_labels = nan(size(datamtx(:,1)));
all_trial_labels(stage_idx==3) = trial_label;
num_trials = length(unique(all_trial_labels(~isnan(all_trial_labels))));
trial_time_bounds = nan(num_trials,2);
for itrial = 1:num_trials 
    trial_time_bounds(itrial,:) = [nanmin(trial_pos(trial_label==itrial,1)) nanmax(trial_pos(trial_label==itrial,1))];
end

%add column
datamtx = [datamtx nan(size(datamtx,1),1)];
for itrial = 1:size(trial_time_bounds,1)
    datamtx(datamtx(:,1)>=trial_time_bounds(itrial,1) & datamtx(:,1)<trial_time_bounds(itrial,2), end) = itrial;
end


%% Limit datamtx to time range

%find frame range
frm_rng = datamtx(:,5);
frm_rng(frm_rng~=0)=nan;
frm_rng(frm_rng==0) = 1:sum(datamtx(:,5)==0);
time_range_idx = datamtx(:,1)>=time_range(1) & datamtx(:,1)<=time_range(2);
frm_rng = frm_rng(time_range_idx);
frm_rng = [nanmin(frm_rng), nanmax(frm_rng)];
datamtx_all = datamtx;
datamtx = datamtx(time_range_idx, :);


%% OUTPUT
datamtx_txyhs = datamtx;


end






function fixed_timestamps = fix_time(timestamps, start_time)
    timestamps = timestamps - start_time;
    fixed_timestamps = timestamps./1000000;
end
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
    
    %index to delete dubious points
    deletions = logical(deletions);
    timeXYD(deletions, 2:3) = NaN;

    %replace deleted points with interpolated values
    non_nan_pos = ~isnan(timeXYD(:,2)) & ~isnan(timeXYD(:,3)); %index
    timeXYD(:,2) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 2), timeXYD(:,1), 'linear');
    timeXYD(:,3) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 3), timeXYD(:,1), 'linear');
    
end
function timeXYD = smooth_pos(timeXYD, smoothing_window_size)
    non_nan_pos = ~isnan(timeXYD(:,2)) & ~isnan(timeXYD(:,3)) & ~isnan(timeXYD(:,4)); %index
    timeXYD(non_nan_pos,2) = smooth(timeXYD(non_nan_pos, 2), smoothing_window_size);
    timeXYD(non_nan_pos,3) = smooth(timeXYD(non_nan_pos, 3), smoothing_window_size);
    timeXYD(non_nan_pos,4) = smooth_hd(timeXYD(non_nan_pos, 4), smoothing_window_size');
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
        display([num2str(sum(EventIDs==4)) ' manual events found. Only using the first 4.'])
    end

    %parse datamtx using 4 manual flags (end bowl, start maze, end maze,
    %start bowl)
    stage_idx = nan(size(datamtx,1), 1);
    manual_flag_times = Timestamps_events(EventIDs==4);
    for ievent = 1:4
       stage_idx(datamtx(:,1)<= manual_flag_times(ievent) & isnan(stage_idx)) = ievent;
       if ievent == 4
           stage_idx(datamtx(:,1) > manual_flag_times(ievent)) = ievent+1;
       end
    end
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
function timeXYD = interp_at_nans(timeXYD)
%replace deleted points with interpolated values
    non_nan_pos = ~isnan(timeXYD(:,2)) & ~isnan(timeXYD(:,3)); %pos index
    non_nan_hd = ~isnan(timeXYD(:,4)); %hd index
    timeXYD(:,2) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 2), timeXYD(:,1), 'linear');
    timeXYD(:,3) = interp1(timeXYD(non_nan_pos, 1), timeXYD(non_nan_pos, 3), timeXYD(:,1), 'linear');
    
    if ~isempty(timeXYD(non_nan_hd, 1))
        timeXYD(:,4) = interp_hd(timeXYD(non_nan_hd, 1), timeXYD(non_nan_hd, 4), timeXYD(:,1), 'linear');
    end
    
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