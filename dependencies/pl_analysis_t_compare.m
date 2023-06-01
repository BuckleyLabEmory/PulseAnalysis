function dt = pl_analysis_t_compare(TimeAxis_DCS,raw_ABP_masterdata,Dbfit,numnans,varargin)
% This function performs the pre-processing first tested in 
% /Users/taraurner/BuckleyLab/CrCp/CrCP_tmu/tmu_testing/MatchABPandDCS_091819_3_tmu.m
%%
    p = inputParser;
    addRequired(p,'TimeAxis_DCS')
    addRequired(p,'raw_ABP_masterdata')
    addRequired(p,'Dbfit')
    addRequired(p,'numnans')
    addOptional(p,'FileID',[])
    parse(p,TimeAxis_DCS,raw_ABP_masterdata,Dbfit,numnans,varargin{:});
    
    
    dcs_data = Dbfit;

    %%
    dt = struct();

    %% Make sure this data has all unique timepoints (not concatenated acquisitions)

    % DCS

    [t_DCS, t_DCS_u, ~] = unique(TimeAxis_DCS);
    
    if length(t_DCS)~=length(t_DCS_u)
        error('Non-unique DCS timepoints')
    end

    [t_ABP, t_ABP_u, ~] = unique(raw_ABP_masterdata(:,1));
    
    % add to output structure
    dt.t_DCS = t_DCS;
    dt.t_ABP = t_ABP;
    
    % This information is written out to the text file output each time
    % pulse_analysis_main is run
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['\nDbfit length = ' num2str(length(Dbfit))...
            '\nTimeAxis_DCS length = ' num2str(length(TimeAxis_DCS))...
            '\nUnique TimeAxis_DCS length = ' num2str(length(TimeAxis_DCS))]);
    else
        fprintf(['\nDbfit length = ' num2str(length(Dbfit))...
            '\nTimeAxis_DCS length = ' num2str(length(TimeAxis_DCS))...
            '\nUnique TimeAxis_DCS length = ' num2str(length(TimeAxis_DCS))])        
    end

    %% Sampling Constants
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['\n\nCalculating sampling constants...\n\n']);
    else    
        fprintf(['\n\nCalculating sampling constants...\n\n']);
    end
    % DCS

    DCS_BFI = dcs_data; 


    dT_all_DCS = diff(t_DCS);

    dT_DCS = mean(dT_all_DCS);
    f_DCS = 1/dT_DCS;
    
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['DCS sampling frequency is: ' num2str(f_DCS) ' hz. Or every '...
        num2str(dT_DCS) ' seconds\n']);
    else
        fprintf(['DCS sampling frequency is: ' num2str(f_DCS) ' hz. Or every '...
        num2str(dT_DCS) ' seconds\n']);
    end
    
    % add to output structure
    dt.f_DCS = f_DCS;
    
    % ABP (A-Line)

    temp = raw_ABP_masterdata(:,3);
    ABP_mmHg = temp(t_ABP_u);

    dT_all_ABP = diff(raw_ABP_masterdata(:,1));
    dT_ABP = mean(dT_all_ABP);
    f_ABP = 1/dT_ABP;


    % add to output structure 
    dt.f_ABP = f_ABP;
    
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['ABP sampling frequency is: ' num2str(f_ABP) ' hz. Or every '...
        num2str(dT_ABP) ' seconds\n']);
    else
        fprintf(['ABP sampling frequency is: ' num2str(f_ABP) ' hz. Or every '...
        num2str(dT_ABP) ' seconds\n']);     
    end


    %% Overlapping signals in time
    
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['\n\nTemporally overlapping signals...\n\n']);     
    else
        fprintf(['\n\nTemporally overlapping signals...\n\n']);
    end
    % ABP signal temporal params

    % check for any NaNs at the beginning or end and drop them
    abp_vals_idx = find(~isnan(ABP_mmHg));
    start_ABP = t_ABP(abp_vals_idx(1));
    end_ABP = t_ABP(abp_vals_idx(end));
    t_ABP_corr = t_ABP(abp_vals_idx(1):abp_vals_idx(end));
    
    length_t_ABP = end_ABP - start_ABP;
        
    
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['The ABP signal starts at ' num2str(start_ABP) ...
            ' seconds and ends at ' num2str(end_ABP) ' seconds'...
            ' for a total length of : ' num2str(length_t_ABP) ' seconds'...
            ' or ' num2str(length_t_ABP/60) ' minutes\n']);         
    else
        fprintf(['The ABP signal starts at ' num2str(start_ABP) ...
            ' seconds and ends at ' num2str(end_ABP) ' seconds'...
            ' for a total length of : ' num2str(length_t_ABP) ' seconds'...
            ' or ' num2str(length_t_ABP/60) ' minutes\n']);
    end
    
    % DCS signal temporal params
    
    % check for any NaNs at the beginning or end and drop them
    dcs_vals_idx = find(~isnan(DCS_BFI));
    start_DCS = t_DCS(dcs_vals_idx(1));
    end_DCS = t_DCS(dcs_vals_idx(end));
    t_DCS_corr = t_DCS(dcs_vals_idx(1):dcs_vals_idx(end));
    DCS_BFI_corr = DCS_BFI(dcs_vals_idx(1):dcs_vals_idx(end));
    
    length_t_DCS = end_DCS - start_DCS;
    
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['The DCS signal starts at ' num2str(start_DCS) ...
            ' seconds and ends at ' num2str(end_DCS) ' seconds'...
            ' for a total length of : ' num2str(length_t_DCS) ' seconds'...
            ' or ' num2str(length_t_DCS/60) ' minutes\n']);        
    else
    fprintf(['The DCS signal starts at ' num2str(start_DCS) ...
        ' seconds and ends at ' num2str(end_DCS) ' seconds'...
        ' for a total length of : ' num2str(length_t_DCS) ' seconds'...
        ' or ' num2str(length_t_DCS/60) ' minutes\n']);
    end
    

    dt.corr_ABP_idx = [abp_vals_idx(1) abp_vals_idx(end)];
    dt.corr_DCS_idx = [dcs_vals_idx(1) dcs_vals_idx(end)];
    dt.corr_ABP_bnds = [start_ABP end_ABP];
    dt.corr_DCS_bnds = [start_DCS end_DCS];
    dt.t_DCS_corr = t_DCS_corr;
    dt.t_ABP_corr = t_ABP_corr;
    dt.DCS_corr = DCS_BFI_corr;
    dt.ABP_corr = ABP_mmHg(abp_vals_idx(1):abp_vals_idx(end));
    
 
    % Which signal starts first in time

    first_signal = min([start_DCS, start_ABP]);

    % If ABP measurement start time leads DCS measurement start time
    if first_signal == start_ABP
        starter_signal = 'ABP';
        delta_ABP_DCS = start_DCS - start_ABP;
        if any(p.Results.FileID)
            fprintf(p.Results.FileID,['The ABP signal starts first by ' num2str(delta_ABP_DCS) '\n']);    
        else
            fprintf(['The ABP signal starts first by ' num2str(delta_ABP_DCS) '\n'])
        end
        dt.delta_ABP_DCS = delta_ABP_DCS;
        dt.delta_DCS_ABP = [];
    % If DCS measurement start time leads ABP measurement start time
    elseif first_signal == start_DCS
        starter_signal = 'DCS';
        delta_DCS_ABP = start_ABP - start_DCS;
        dt.delta_DCS_ABP = delta_DCS_ABP;
        dt.delta_ABP_DCS = [];
        if any(p.Results.FileID)
            fprintf(p.Results.FileID,['The DCS signal starts first by ' num2str(delta_DCS_ABP) '\n']);
        else
            fprintf(['The DCS signal starts first by ' num2str(delta_DCS_ABP) '\n'])
        end    
    end

    dt.first_signal = first_signal;
    dt.starter_signal = starter_signal;

    
    %% Resample ABP at BFI timepoints


    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['\n\nDownsampling to a constant frequency...\n\n']);        
    else 
        fprintf(['\n\nDownsampling to a constant frequency...\n\n']);
    end

    t_ds_start = dcs_vals_idx(1);
    t_ds_end = dcs_vals_idx(end);
    dt_ds = dT_DCS;
    t_ds = t_DCS_corr;
    sampling_f = f_DCS;

    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['Downsampling time interval is from ' num2str(t_ds_start) ...
        ' to ' num2str(t_ds_end) ' in ' num2str(dt_ds) ' second intervals.\n'])
    else
        fprintf(['Downsampling time interval is from ' num2str(t_ds_start) ...
        ' to ' num2str(t_ds_end) ' in ' num2str(dt_ds) ' second intervals.\n'])
    end
    
    % Calculate downsample rate ABP

    ABP_downsample_f = round(f_ABP/sampling_f);

    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['Downsample ABP signal at a rate of approximately '...
        num2str(ABP_downsample_f) ' hz for a sampling rate of ' ...
        num2str(sampling_f) ' hz \n'])        
        fprintf(p.Results.FileID,['Downsampling time interval is from ' num2str(t_ds_start) ...
        ' to ' num2str(t_ds_end) ' in ' num2str(dt_ds) ' second intervals.\n'])
    else
        fprintf(['Downsample ABP signal at a rate of approximately '...
        num2str(ABP_downsample_f) ' hz for a sampling rate of ' ...
        num2str(sampling_f) ' hz \n'])        
        fprintf(['Downsampling time interval is from ' num2str(t_ds_start) ...
        ' to ' num2str(t_ds_end) ' in ' num2str(dt_ds) ' second intervals.\n'])
    end

    % Interpolate ABP

    ABP_downsample = interp1(t_ABP_corr,ABP_mmHg(abp_vals_idx(1):abp_vals_idx(end)),...
        t_ds);

    dt.ABP_downsample = ABP_downsample;
    
    if any(p.Results.FileID)
        fprintf(p.Results.FileID,['Both signals now sampled from t = ' num2str(t_ds(1)) ' to ' ...
        num2str(t_ds(end))])    
    else
        fprintf(['Both signals now sampled from t = ' num2str(t_ds(1)) ' to ' ...
        num2str(t_ds(end))])        
    end
    
 %% Interp BFI points below num nan thresh
 
    DCS_vals = ~isnan(DCS_BFI_corr);
    DCS_nans = isnan(DCS_BFI_corr);  
    [B,N,bIB] = RunLength(DCS_BFI_corr);
    over_thresh = N > numnans;
    
    
    idx_over_thresh = bIB(over_thresh);
    n_over_thresh = N(over_thresh) - 1;
    interp_t = ones(length(DCS_BFI_corr),1);
    if ~isempty(idx_over_thresh)
       interp_t(idx_over_thresh:idx_over_thresh+n_over_thresh) = 0;
    
    end
    interp_t = logical(interp_t);
    interp_nans = DCS_nans & interp_t;
    BFI_interp_vals = interp1(t_DCS_corr(DCS_vals),DCS_BFI_corr(DCS_vals),t_DCS_corr(interp_nans));
   
    BFI_interp_nativef = DCS_BFI_corr;
    BFI_interp_nativef(interp_nans) = BFI_interp_vals;
    dt.BFI_interp_nativef = BFI_interp_nativef;    
 
 
 %% Resample BFI at ABP timepoints
    
    abp_overlap = (dt.t_ABP_corr >= dt.t_DCS_corr(1) & dt.t_ABP_corr <= dt.t_DCS_corr(end));
    dt.ABP_overlap = abp_overlap;
    BFI_interp = interp1(dt.t_DCS_corr,dt.BFI_interp_nativef,dt.t_ABP_corr(abp_overlap),'pchip');
    dt.BFI_interp = BFI_interp;
    dt.t_BFI_interp = dt.t_ABP_corr(abp_overlap);
    dt.ABP_binterp = dt.ABP_corr(abp_overlap);
    
end