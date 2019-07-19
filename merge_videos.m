%merge video files

%% load

%filepath
fp = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\videos\m1_6-15_trial-';
%file names
fnames = {'05', '06', '07', '08', '09', '10', '11', '12', '13'};
%filetype
ft = '.avi';


%load and merge video files
vfs = [];
for i = 1:length(fnames)
    
    vfs = [vfs; load_video_frames([fp fnames{i} ft])];    
end


%% save

%save details for new video
save_filename = 'merged_video';
clock_hold = clock;
clock_hold = [num2str(clock_hold(1)) '-' num2str(clock_hold(2)) '-' num2str(clock_hold(3))];
outfile_vid = ['C:\Users\ampm1\Documents\MATLAB\tt_ephys\videos\' save_filename '_' clock_hold '.avi'];

%preallocate new video
video_out = VideoWriter(outfile_vid);
v.Quality = 95;
video_out.FrameRate = (29)*10;


%iteratively edit each vid_frame and save a video
open(video_out)
figure;
for ivf = 1:length(vfs)
    
    %plot frame
    imagesc(vfs{ivf})
    
    %aesthetics
    set(gca,'TickLength',[0, 0])
    axis equal
    axis tight
    axis off
    set(gcf, 'Position', [813,589,1125,749])
      
    %create image from plot
    frame = getframe;

    %write video frame from image
    writeVideo(video_out, frame.cdata);

end
close(video_out)