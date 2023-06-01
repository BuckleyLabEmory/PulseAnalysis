function ftbl =pl_analysis_ftbl_gen(patient,state,t_abp_raw,abp_raw,fs_abp,abp_downsample,t_bfi,bfi,fs_bfi,t_bfi_interp,bfi_interp,abp_binterp)
% This code just adds all of the time series data to one row for each subject in one big table, 
% making it easier to pass this into the next steps
% 't_abp_raw','abp_raw','fs_abp','t_bfi','bfi','fs_bfi' --> raw time axes and data for ABP and BFI
% 'abp_downsample', --> ABP interpolated to BFI time axis
% 't_bfi_interp', --> 125 ABP time axis
% 'bfi_interp', --> BFI interpolated to ABP time axis
% 'abp_binterp'} --> this is exactly the same as raw abp
    ftbl = [patient,state,{t_abp_raw},{abp_raw},fs_abp,{t_bfi},{bfi},fs_bfi,{abp_downsample},{t_bfi_interp},{bfi_interp},{abp_binterp}];
    ftbl = cell2table(ftbl);
    ftbl.Properties.VariableNames = {'name','state',...
        't_abp_raw','abp_raw','fs_abp','t_bfi','bfi','fs_bfi','abp_downsample','t_bfi_interp','bfi_interp','abp_binterp'};
end