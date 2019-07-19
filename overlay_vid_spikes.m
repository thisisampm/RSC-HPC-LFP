%function video_out = overlay_vid_spikes(folderpath, time_range, clusters, savename, ripple_idx, csc_mtx)
%plot spikes overlaying raw nlx video
%{
%load position and spike data for time range
[datamtx_tr, frm_rng, datamtx_all] = loadnlx_raw(folderpath, time_range);
datamtx_vidonly = datamtx_tr(datamtx_tr(:,5)==0, :);

%hard correct weird nlx error
frm_rng = frm_rng - 5;

%load video frames
vid_frames = load_video_frames([folderpath '\VT1.mpg'], frm_rng);

%get matlab colors
figure;
color_mtx = get(gca,'ColorOrder');
close
%}
%save details for new video
save_filename = savename;
outfile_vid = ['C:\Users\ampm1\Documents\MATLAB\tt_ephys\videos\' save_filename '.avi'];

%preallocate new video
video_out = VideoWriter(outfile_vid);
v.Quality = 100;
video_out.FrameRate = sum(datamtx_tr(:,5)==0)/...
    (max(datamtx_tr(datamtx_tr(:,5)==0,1)) - min(datamtx_tr(datamtx_tr(:,5)==0,1)));
video_out.FrameRate = (video_out.FrameRate)/6;

%iteratively edit each vid_frame and save a video
open(video_out)
figure;
for ivf = 1:length(vid_frames)
    
    %plot frame
    imagesc(vid_frames{ivf})
    
    %aesthetics
    set(gca,'TickLength',[0, 0])
    axis equal
    axis tight
    axis off
    %set(gcf, 'Position', [813,589,1125,749])
    set(gcf, 'Position', [129         374        1111         604])
    
    
    %ADD EDITS
    %
    hold on

        %plot all spikes occuring before this frame
        for iclust = 1:size(clusters,1)
            %current cluster
            clust = clusters(iclust,1);

            %cluster (spike event) index
            clust_idx = datamtx_all(:,5)==clust;

            %time index & TRIAL INDEX
            spike_time_idx = datamtx_all(:,1) <= datamtx_vidonly(ivf,1) & datamtx_all(:,6)>0 & datamtx_all(:,6)<=13;

            %plot
            plot(datamtx_all(clust_idx & spike_time_idx, 2), datamtx_all(clust_idx  & spike_time_idx, 3), 'g.', 'markersize', 20)%, 'color', color_mtx(1,:))
        end

        %plot current position
        spike_time_idx = datamtx_tr(:,1) <= datamtx_vidonly(ivf,1);
        plot(datamtx_vidonly(ivf,2), datamtx_vidonly(ivf,3), 'y+', 'markersize', 20)
        
        
        %plot lfp
        for col = 2

            %time index
            datamtx_vidonly(ivf,1)
            
            lo_time = datamtx_vidonly(ivf,1) - .5;
            hi_time = datamtx_vidonly(ivf,1);
            lfp_time_idx = datamtx_all(:,1) >= lo_time & datamtx_all(:,1) <= hi_time;
            
            %position
            vidwidth = 100;
            vidheight = 100;
            xpos = vidwidth*(2/3) + 270;
            ypos = vidheight*(1/3)-115;
            
            %size
            xpos_size = vidwidth*(1/3) + 250;
            ypos_size = vidheight*(1/6) - 5;
            
            %xaxis
            cscidx = csc_dmidx(csc_mtx, datamtx_all, lfp_time_idx);
            xplot = norm_mtx((1:sum(cscidx))').*xpos_size + xpos;
            
            %yaxis
            norm_csc = csc_mtx(:,col)./10000; %shrink (not norm)
            yplot = (norm_csc(cscidx).*ypos_size) + ypos;
            yplot = -yplot;
            
            %ripple index
            rip_idx = csc_dmidx(csc_mtx, datamtx_sesh, ripple_idx);
            
            xplot_select = find(rip_idx(cscidx)==1);
            xplot_deselect = find(rip_idx(cscidx)~=1);
            %xplot_rip = xplot(xplot_select);
            %yplot_rip = yplot(xplot_select);
            
            xplot_rip = xplot; xplot_rip(xplot_deselect) = nan;
            yplot_rip = yplot; yplot_rip(xplot_deselect) = nan;
           
            %plot
            plot(xplot, yplot, 'b-', 'linewidth', 1)
            plot(xplot_rip, yplot_rip, 'r-', 'linewidth', 1)
 
        end
        
        
        %PLOT SPIKES IN TIME WITH LFP
        lfp_window_duration = hi_time - lo_time;
        clust_idx = datamtx_all(:,5)>=6;
        lfp_spk_times = datamtx_all(lfp_time_idx & clust_idx,1);
        lfp_spk_times = (lfp_spk_times - lo_time)/lfp_window_duration;
        lfp_spk_times = lfp_spk_times.*xpos_size + xpos;

        lfp_spk_ypos = repmat(ypos+220, size(lfp_spk_times));
        plot(lfp_spk_times, lfp_spk_ypos, 'g.', 'markersize', 20)
        

    hold off
      
    %create image from plot
    pause(0.1)
    frame = getframe;

    %write video frame from image
    writeVideo(video_out, frame.cdata);

end
close(video_out)


