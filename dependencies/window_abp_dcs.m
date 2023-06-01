function [window_start, window_end] = window_abp_dcs(t_ds,sampling_f,window_size)
% 
% window_num = floor((pdo.t_ds(end)-pdo.t_ds(1))/window_size);
% 
% fprintf(['\n\nWindowing ' num2str(pdo.t_ds(end)) ' sec long signal in ' ...
%     num2str(window_size) ' increments for a total of ' num2str(window_num) ' windows...\n\n'])
% 
% fcs = (1:window_num).*(window_size+pdo.t_ds(1));
% 
% window_start = zeros(window_num,1);
% window_end = zeros(window_num,1);
% 
% window_start(1) = 1;
% 
% for a=1:window_num
%     idx = find(abs(pdo.t_ds - fcs(a)) == min(abs(pdo.t_ds - fcs(a))));
%     window_end(a) = idx;
%     window_start(a+1) = idx+1;
%     fprintf(['Window length ' num2str(length(window_start(a):window_end(a))) '\n'])
% end
%%
    window_length = ceil(window_size/(1/sampling_f));
    window_num = floor(length(t_ds)/window_length);

    windows = 0:window_length:(window_length*window_num)+1;
    window_start = (windows(1:end-1)+1)';
    window_end = (windows(2:end)+1)';

end