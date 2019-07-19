function [vid_frames, vid_frames_times] = load_video_frames(filepath, varargin)
%load raw nlx video

if nargin == 2
    frame_range = varargin{1};
else
    frame_range = [0 inf];
end


%create vid object
v = VideoReader(filepath);

%preallocate video frame time vector
vid_frames_times = 0;

%createa movie structure array
vidStruct = struct('cdata', zeros(v.Height, v.Width, 3, 'uint8'), 'colormap', []);

%read frames from start time to time_int(2)
k = 1; k2 = 0;
while hasFrame(v)
    if k >= frame_range(1) && k <= frame_range(2)
        k2 = k2+1;
        vidStruct(k2).cdatat = readFrame(v);
    elseif k > frame_range(2)
        break
    else
        readFrame(v);
    end
    k = k+1;
end


%convert to cell array
vid_frames = struct2cell(vidStruct);
vid_frames = permute(vid_frames, [3, 2, 1]); 

%video frames only
vid_frames = vid_frames(:,:,3);

