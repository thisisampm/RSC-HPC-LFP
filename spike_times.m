function [st_out, st_out_norm] = spike_times(datamtx, clusters, time_bounds, varargin)
%st out is a cell length(clusters) containing the times of spikes
%for each cell during time_bounds (inclusive or exclusive)

%exclusive or inclusive
if nargin > 3
    inclusive = varargin{1};
else
    inclusive = 1;
end
%preallocate cell
st_out = cell(size(clusters,1),1);
st_out_norm = cell(size(clusters,1),1);

%for each cell
for iclust = 1:size(clusters,1)
    
    cluster_idx = datamtx(:,8)==clusters(iclust,1);
    if inclusive==1
        time_idx = datamtx(:,1)>=time_bounds(1) & datamtx(:,1)<=time_bounds(2);
    elseif inclusive==0
        time_idx = datamtx(:,1)>time_bounds(1) & datamtx(:,1)<time_bounds(2);
    end
    st_out{iclust} = datamtx(cluster_idx & time_idx, 1);
    st_out_norm{iclust} = (st_out{iclust}-time_bounds(1))/diff(time_bounds);
    
end


