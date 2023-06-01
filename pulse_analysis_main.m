function [ftbl,wtbl,ptbl,dat_tbl] = pulse_analysis_main(data,savedir,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Copyright (c) 2023 Tara Urner              |\__/,|   (`\ 
%  http://buckleylab.gatech.edu/            _.|o o  |_   ) )
%                                          -(((---(((--------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use ABP waveform at native frequency to identify ABP cardiac pulses and
% pull out overlapping BFI data.
% 
% Inputs:
% Required: 
% data: Folder of input data in the expected format.
% savedir: directory to save output text file
% Optional:
% 'sds' = source detector separation to consider. Doesn't change the
%         analysis at all, just a useful column to have in the output
%         table. If none given, empty column is created.
% 'corr_method':
%   'xcorr': Default, shift windows by maximizing correlation returned by
%            xcorr
%   'corrcoef': Other option, shift windows by maximizing correlation
%               returned by corrcoef
% 'WindowType': 
%   'All': Window entire timeseries with no overlaps. This is actually the
%           only fully developed method right now, and so is the default. Left in
%           as an input for future development.
% 'WindowSize': Windows size in seconds, 15s default
% 'FileText': Text to save with the data, in this example the "run" script
%               that calls this function
% 'AnalyzeMarks': Optional table input is given to sort marks by state. See
%                 run_pulse_analysis_usemarks1 or 2.
% 'NumInterpNaNs': Number of of consecutive NaNs in BFI data allowed before
%                  window is dicarded. Default 2.
% 'cnap': if the CNAP was used for acquisition.

%% Parse inputs

    %%% DEFAULTS %%%
    defaultWindowAll = 'All';
    defaultWindow = 15; % seconds
    defaultNaNs = 2;
    %defaultInterp = 100;
    defaultCorr = 'xcorr';
    %%%%%%%%%%%%%%%%
    
    p = inputParser;
    addRequired(p,'data');
    addRequired(p,'savedir');
    addOptional(p,'sds',[]);
    addOptional(p,'corr_method',defaultCorr);
    addOptional(p,'WindowType',defaultWindowAll);
    addOptional(p,'WindowSize',defaultWindow);
    addOptional(p,'FileText',[]);
    addOptional(p,'AnalyzeMarks',[]);
    addOptional(p,'NumInterpNaNs',defaultNaNs)
    addOptional(p,'cnap',[])
    
    parse(p,data,savedir,varargin{:});

    fprintf(['RUNNING WITH ' p.Results.corr_method ' CORRELATION METHOD'])

    
%% Run table

% load analysis table
%run_tbl_dir = p.Results.AnalyzeMarks;
% [run_tbl,~] = dcs_dir_load(run_tbl_dir);
% run_sess = categorical(run_tbl.SubjectID); 
    
    %% Analysis record

psum = struct();
asum = struct();

%% Parse DB
asum.inputs = p.Results;
asum.db_dir = p.Results.data;


[asum.load_table,asum.patient_db] = parse_db(asum.db_dir);


if ~isempty(p.Results.AnalyzeMarks)
        run_tbl_dir = p.Results.AnalyzeMarks;
        run_tbl = readtable(run_tbl_dir);
        run = run_tbl.Run;
        %[run_tbl,run] = dcs_dir_load(run_tbl_dir);
        run_tbl = run_tbl(logical(run),:);
        run_sess = categorical(run_tbl.SubjectID);
        select = ismember(categorical(asum.patient_db(:,1)),run_sess);
        asum.patient_db = asum.patient_db(select,:);
end
%% Save Dir

asum.save_dir_base = p.Results.savedir;

t = datetime('now');
t.Format = 'yy-MM-dd_HH-mm-ss';

outfile = [asum.save_dir_base char(t) '_run_out.txt'];
fileID = fopen(outfile,'w');

if any(p.Results.FileText)
    fprintf(fileID,p.Results.FileText)
    fprintf(p.Results.FileText,['\n\nParsing patient database from ' asum.db_dir '...\n\n'])
else
    fprintf(['\n\nParsing patient database from ' asum.db_dir '...\n\n'])
end

%% Parse Patients
patients = [char(asum.patient_db{:,1})];
%sessions = cellfun(@(x) size(x,1),{asum.patient_db{:,3}});
for i=1:size(patients,1)
  psum(i).patient = string(patients(i,:));
%  psum(i).num_sessions = sessions(i);
end 

num_p = length(psum);
%total_sess = sum([psum.num_sessions]);

if any(p.Results.FileText)
    fprintf(fileID,p.Results.FileText)
end

fprintf(fileID,['\n\nParsing patient database from ' asum.db_dir '...\n\n'])

% Table and struct init

pdo = struct();
ftbl = table();
wtbl = table();
ptbl = table();
dat_tbl = table();
wv = table();
rc = 1;

%% The big processing loop for every subject
for i=1:num_p
    patient_i = asum.patient_db{i,1} % which subject is this
    pdo(rc).patient = patient_i;

    if ~isempty(p.Results.AnalyzeMarks)
        run_tbl_tmp = run_tbl(run_sess == patient_i,:);
        if isempty(run_tbl_tmp) % If the code can't find data for this subject in the input folder
            error(['No data matching subject id ' patient_i ' present in ' data])
        end
        run_tbl_in = run_tbl_tmp;
        state = 'Defined' % to be determined per window/pulse based on marks table
    else
        state = 'Not Defined' % no marks table provided
    end  

    % load the data
    numvar = size(asum.patient_db{i,3}{3},1);
    for k=1:numvar
        load(asum.patient_db{i,3}{3}.path(k));
    end 
    %end
    
    if ~isempty(p.Results.cnap)
       % Converts voltage to ABP for the CNAP
       raw_ABP_masterdata(:,3) = ((raw_ABP_masterdata(:,2)./5)*500);
    end    

    if ~isempty(p.Results.sds) & size(Dbfit,2) > 1
    % if the sds argument is used and Dbfit is a column vector with
    % multiple columns, assumes sds is the column number to use.
    Dbfit = Dbfit(:,p.Results.sds);
    end
    
    %% Resampling, interpolating single NaNs by default
    numnans = p.Results.NumInterpNaNs;
    dt = pl_analysis_t_compare(TimeAxis_DCS,raw_ABP_masterdata,...
        Dbfit,numnans,'FileID',fileID);  
    
    %% Full time series table
    ftbl_n =pl_analysis_ftbl_gen(patient_i,state,dt.t_ABP_corr,dt.ABP_corr,dt.f_ABP,dt.ABP_downsample,...
        dt.t_DCS_corr,dt.DCS_corr,dt.f_DCS,dt.t_BFI_interp,dt.BFI_interp,dt.ABP_binterp);
    
    ftbl_n.dt = dt; % sampling information for raw ABP and BFI data
    ftbl_n.data_dir = data; % provided path to the data
    ftbl_n.corr_method = {p.Results.corr_method}; % xcorr by default
    if ~isempty(p.Results.AnalyzeMarks)
        ftbl_n.analysis_tbl_dir = {run_tbl_dir};
        ftbl_n.marks = {Marks_index};
        ftbl_n.run_tbl = {run_tbl_in};
    end
    if ~isempty(p.Results.sds)  % if sds was given, add to this table 
        sds = p.Results.sds;
        sds = repmat({sds},height(ftbl_n),1);
        tmp = cell2table(sds);
        ftbl_n = [tmp, ftbl_n];
    else    % Otherwise, sds is an empty column
        sds = repmat({''},height(ftbl_n),1);    
        ftbl_n.sds = sds;
    end
    ftbl = [ftbl ; ftbl_n];
    
    %% Create windows
    window_size = p.Results.WindowSize; % 15s by default
    
    %%% WINDOWING %%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Window entire signal with no overlap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(p.Results.WindowType,'All') % written to support adding other types of windowing
        
             [window_start, window_end] ...
             = window_abp_dcs(ftbl_n.t_bfi{:},ftbl_n.fs_bfi,...
             window_size);    
    end
    
    %% Sort windows into states - either defined by marks or 'none'.
    wtbl_n = window_state_sort(ftbl_n,window_start,window_end);

    %%%%%%%%%%%%%%% Most generalized pulse analysis %%%%%%%%%%%%%%%%%%%%%%%
    % For anyone wanting to insert this pulse analysis method 
    % into an existing pipeline that handles all of the above 
    % loading/syncing/windowing for a different experimental setup, 
    % these three subfunctions are what you want to use.
    %% Extract the window data, identify ABP pulses, calculate shift between ABP and BFI
    wtbl_n = pl_analysis_interp_wtbl_gen_v2(ftbl_n,wtbl_n,numnans);

    %% Extract BFI pulses based on ABP reference, add to pulse table
    pl_n = pl_analysis_ptbl_gen_v3(ftbl_n,wtbl_n);
    ptbl = [ptbl ; pl_n];

    %% Generate window-averaged pulses and analyze pulse morphology
    % If you want to change bin size that is hard-coded in here
    [dat_tbl_n,wtbl_n] = dat_tbl_gen_v3(wtbl_n,pl_n);
    wtbl = [wtbl ; wtbl_n];
    dat_tbl = [dat_tbl ; dat_tbl_n];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    clear TimeAxis_DCSraw_ABP_masterdata Marks_index Dbfit    
end

end