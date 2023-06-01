function peaks = getpeaks3(abp,sf)   
% This code calls the delineator function to identify pulses
    % Reference:
    %   BN Li, MC Dong & MI Vai (2010) 
    %   On an automatic delineator for arterial blood pressure waveforms
    %   Biomedical Signal Processing and Control 5(1) 76-81.
    
    % LI Bing Nan @ University of Macau, Feb 2007
    %   Revision 2.0.5, Apr 2009
% Various indices are then computed and stored in the 'peaks' structure

    peaks = struct;

    [onsetp,peakp,dicron] = delineator(abp,sf);
    
    % raw
    peaks.ao_r = onsetp;
    peaks.ap_r = peakp;
    peaks.ad_r = dicron;

    eidx = min([length(peaks.ao_r) length(peaks.ap_r) length(peaks.ad_r)]);
    
    peaks.ao = peaks.ao_r(1:eidx); % pulse onset
    peaks.ap = peaks.ap_r(1:eidx); % pulse peak
    peaks.ad = peaks.ad_r(1:eidx); % pulse dicrotic notch
       
    
    bp_pulse = [peaks.ao;peaks.ap;peaks.ad]; % indices of all features in one array
    peaks.a_p = bp_pulse;    


    % X-axis (time) features
    t = struct;
    % pressure
    t.aoo = diff(peaks.ao); % diff between onsets
    t.aop = peaks.ap - peaks.ao; % diff between onsets and peaks
    t.aod = peaks.ad - peaks.ao; % diff between dicrotic notch and onsets
    t.app = diff(peaks.ap); % peak-to-peak diff
    t.apd = peaks.ad - peaks.ap;    % dicrotic notch to peak diff
    t.add = diff(peaks.ad); % dicrotic notch-to-notch diff
    
    peaks.t = t; % store in t field - for time axis features
    
    % Y axis (pulse units) features
    y = struct;
    y.aop = abp(peaks.ap) - abp(peaks.ao); % pulse onset to peak height
    y.aod = abp(peaks.ad) - abp(peaks.ao); % pulse onset to dicrotic notch height
    y.adp = abp(peaks.ap) - abp(peaks.ad); % pulse dicrotic notch to peak height
    y.app = diff(abp(peaks.ap)); % difference in peak values
    y.add = diff(abp(peaks.ad)); % difference in dicrotic notch values
    y.aoo = diff(abp(peaks.ao)); % difference in pulse onset values

    peaks.y = y; % store in y field - for y axis features

end
