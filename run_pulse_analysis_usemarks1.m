%%%%%%%%%%%%%%%%%%
% Set filepaths  %                                                         
%%%%%%%%%%%%%%%%%%
BASEPATH =  'Y:\buckley-lab\Projects\Waveform_analysis\0_Papers\2023_PulsatilityPaper\000_BIN\FOR GITHUB\BOE_GithubExample'
data_folder = [BASEPATH filesep 'ExampleData'];   % BFI and ABP data to be analyzed
analysis_tbl = [BASEPATH filesep 'data_index_usemarks1.xlsx']  % Contains marks input
save_dir = [BASEPATH filesep 'ExampleOutput' filesep 'marks1_pulse_analysis_out' filesep];      % Save location      % Save location
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

%%%%%%%%%%%%%%%%
% Run analysis %
%%%%%%%%%%%%%%%%
sds = 2;
fname = mfilename('fullpath');
ftext = fileread([fname '.m']);
text = [fname ':' ftext];
[ftbl,wtbl,ptbl,dat_tbl] = pulse_analysis_main(data_folder,save_dir,'sds',sds,...
    'FileText',text,'WindowType','All','WindowSize',15,'AnalyzeMarks',analysis_tbl,'corr_method','xcorr','cnap','true');

%%%%%%%%%%%%%
% Save data %
%%%%%%%%%%%%%
%%% dat_tbl is the recommended variable for export
cd(save_dir);
save([save_dir filesep 'dat_tbl.mat'], 'dat_tbl','-v7.3');

%save([save_dir 'full_timseries_table.mat'], 'ftbl','-v7.3');
%save([save_dir 'window_tbl.mat'], 'wtbl','-v7.3'); % I
%save([save_dir 'pulse_table.mat'], 'pl','-v7.3');

%%%%%%%%%%%%%%%%%%%%
% Visualize pulses %
%%%%%%%%%%%%%%%%%%%%
visualize_pulses