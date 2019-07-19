function plotlfptrace(csc_mtx, timerng, columns)
figure; hold on; 

time_idx = csc_mtx(:,1)>timerng(1) & csc_mtx(:,1)<timerng(2);
time = csc_mtx(time_idx,1); 

colmtx = nan(length(columns), length(time));
for i = 1:length(columns)
    colmtx(i,:) = csc_mtx(time_idx,columns(i));
    %colmtx(i,:) = smooth(colmtx(i,:),50);
    plot(time, colmtx(i,:));
end

end