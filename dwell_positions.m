function [dt_out] = dwell_positions(datamtx, position_bounds, varargin)
%st out is a cell length(clusters) containing the times of spikes
%for each cell during time_bounds (inclusive or exclusive)

%exclusive or inclusive
if nargin > 3
    inclusive = varargin{1};
else
    inclusive = 1;
end

%compute ouput
vid_idx = datamtx(:,8)==0;
if inclusive==1
    pos_idx = datamtx(:,2)>=position_bounds(1) & datamtx(:,2)<=position_bounds(2);
elseif inclusive==0
    pos_idx = datamtx(:,2)>position_bounds(1) & datamtx(:,2)<position_bounds(2);
end
dt_out{1} = datamtx(vid_idx & pos_idx, 2);
    



