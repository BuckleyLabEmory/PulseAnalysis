function DCS_resample = interp_nans(DCS_BFI_corr,t_DCS_corr,test)

% Interpolate DCS, interpolate out NaNs -- As a V2, don't interpolate
    % over multiple NaNs.
    
    DCS_vals = ~isnan(DCS_BFI_corr);
    DCS_nans = isnan(DCS_BFI_corr);

    DCS_NaNvals = interp1(t_DCS_corr(DCS_vals), DCS_BFI_corr(DCS_vals),...
       t_DCS_corr(DCS_nans));
   DCS_resample = DCS_BFI_corr;
   DCS_resample(DCS_nans) = DCS_NaNvals;
   
   if test == 1
       test_nan_interp(t_DCS_corr,DCS_resample,DCS_nans)
   end
end