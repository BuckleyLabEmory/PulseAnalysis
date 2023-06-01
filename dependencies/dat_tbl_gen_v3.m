function [dat_tbl,wtbl ] = dat_tbl_gen_v3(wtbl,pl)


wstates = categorical(wtbl.state);
wsds = categorical(wtbl.state);

placehold = cell(height(wtbl),1);
clrs = lines(4);

%%%%%%%%%%%%%%%%%
bin_size = 0.05;
%%%%%%%%%%%%%%%%%

% Allocating memory
% Pulses from each window will be stacked, preserving relative shifts
wtbl.bfi_pstack = placehold;
wtbl.bfi_tstack = placehold;
wtbl.abp_pstack = placehold;
wtbl.abp_tstack = placehold;
wtbl.abp_dicron = placehold;
%wtbl.clr = placehold;
wtbl.abp_peaks = placehold;
wtbl.abp_pk_next = placehold;
wtbl.abp_tpk_next = placehold;

%% PULSE ALIGNMENT
w_empty = table();
empty_idx = zeros(height(wtbl),1);
for w=1:height(wtbl)
        
       wstart = wtbl.wstart(w);
       wshift = wtbl.tshift_maxcorr(w);
       wstate = wtbl.state{w};
       wsds = wtbl.sds(w);
       
       wt_125hz = wtbl.wt_125hz{w};
       
       % if wsds==1 & strfind(wstate,states{1})
       %     wtbl.clr(w) = {clrs(1,:)};
       % elseif wsds==2 & strfind(wstate,states{1})
       %      wtbl.clr(w) = {clrs(2,:)};
       % elseif wsds==1 & strfind(wstate,states{2})
       %     wtbl.clr(w) = {clrs(3,:)};
       % elseif wsds==2 & strfind(wstate,states{2})
       %     wtbl.clr(w) = {clrs(4,:)};
       % end
       
       pl_w = pl(pl.wstart == wstart,:);
       len_abp = cellfun(@(x) length(x),pl_w.abp_pt);
       mlen_abp = max(len_abp);
       len_bfi = cellfun(@(x) length(x),pl_w.bfi_pt);
       mlen_bfi = max(len_bfi);
       bfi_tstack = nan(height(pl_w),mlen_bfi);
       bfi_pstack = nan(height(pl_w),mlen_bfi);
       abp_tstack = nan(height(pl_w),mlen_abp);
       abp_pstack = nan(height(pl_w),mlen_abp);
       abp_dicron = nan(height(pl_w),1);
       abp_peaks = cell(height(pl_w),1);
       abp_pk_next = nan(height(pl_w),1);
       abp_tpk_next = nan(height(pl_w),1);
       
       for p=1:height(pl_w)
           % BFI shifted relative to ABP pulse
            bfi_tstack(p,1:len_bfi(p)) = pl_w.bfi_pt{p} - pl_w.abp_pt{p}(1) + wshift; % bfi shifted by xcorr amount, normalized to ABP onset
            bfi_pstack(p,1:len_bfi(p)) = pl_w.bfi_p{p};
            abp_tstack(p,1:len_abp(p)) = pl_w.abp_pt{p} - pl_w.abp_pt{p}(1); % abp normalized to start at 0
            abp_pstack(p,1:len_abp(p)) = pl_w.abp_p{p};
            abp_dicron(p) = wt_125hz(pl_w.a_dicron(p)) - pl_w.abp_pt{p}(1);
            % Here we find the start of diastolic runoff, defined as the
            % next peak after the dicrotic notch (previously identified by
            % delineator.m)
            [~,pks] = findpeaks(abp_pstack(p,1:len_abp(p)),'MinPeakDistance',0.2);
            tpks = abp_tstack(p,pks);
            abp_peaks(p) = {tpks};
            next_pks = find(tpks > abp_dicron(p)); % next peak after dicrotic notch
            if isempty(next_pks)
               next_pk = NaN;
               next_tpk = NaN;
            else
                next_pk = next_pks(1);
                next_tpk = tpks(next_pk);                
            end
            abp_pk_next(p) = next_pk;
            abp_tpk_next(p) = next_tpk;
       end
       if isempty(bfi_tstack)
           w_empty = [w_empty;wtbl(w,:)];
           empty_idx(w)=1;
           %fprintf('empty')
       end
       wtbl.bfi_tstack(w) = {bfi_tstack};
       wtbl.bfi_pstack(w) = {bfi_pstack};
       wtbl.abp_tstack(w) = {abp_tstack};
       wtbl.abp_pstack(w) = {abp_pstack};
       wtbl.abp_dicron(w) = {abp_dicron};
       wtbl.abp_peaks(w) = {abp_peaks};
       wtbl.abp_pk_next(w) = {abp_pk_next};
       wtbl.abp_tpk_next(w) = {abp_tpk_next};
       
end

%% PULSE NORMALIZATION
wtbl = wtbl(~empty_idx,:);

% Prior to averaging: normalize ABP on [0 1], scale BFI by the same factor
placehold = cell(height(wtbl),1);
wtbl.abp_tnorm = placehold;
wtbl.bfi_tnorm = placehold;
wtbl.abp_dnorm = placehold;
wtbl.abp_pnorm = placehold;
wtbl.abp_avg_l = placehold;
for w=1:height(wtbl)
    abp_tstack = wtbl.abp_tstack{w};
    bfi_tstack = wtbl.bfi_tstack{w};
    abp_pstack = wtbl.abp_pstack{w};
    bfi_pstack = wtbl.bfi_pstack{w};
    abp_dicron = wtbl.abp_dicron{w};
    abp_pnext = wtbl.abp_tpk_next{w};
    % Scaling factor to normalize all pulses prior to averaging
    abp_ends = max(abp_tstack');
    norm_factor = (1./abp_ends)';
    abp_tnorm = abp_tstack.*norm_factor;
    bfi_tnorm = bfi_tstack.*norm_factor;
    abp_dnorm = abp_dicron.*norm_factor; % new position of dicrotic notch
    abp_pnorm = abp_pnext.*norm_factor; % new position of start of diastolic runoff
    wtbl.abp_tnorm(w) = {abp_tnorm};
    wtbl.bfi_tnorm(w) = {bfi_tnorm};
    wtbl.abp_dnorm(w) = {abp_dnorm};
    wtbl.abp_pnorm(w) = {abp_pnorm};
    wtbl.abp_avg_l(w) = {mean(abp_ends)};
%     if w==1
%         figure
%         for p=1:size(abp_tstack,1)
%             bfi_norm = normalize(bfi_pstack(p,:),'range',[0 1]);
%             abp_norm = normalize(abp_pstack(p,:),'range',[0 1]);
%             plot(abp_tnorm(p,:),abp_norm,'.-','Color',[1.00,0.50,0.50],'MarkerSize',10,'LineWidth',0.5)
%             hold on
%             plot(bfi_tnorm(p,:),bfi_norm,'.-','color',[0.70,0.90,1.00],'MarkerSize',10,'LineWidth',0.5)
%             plot(wtbl.abp_dnorm{w}(p),abp_norm(abp_tnorm(p,:) == wtbl.abp_dnorm{w}(p)),'o','MarkerSize',10,'color','black','Linewidth',2)
%             plot(wtbl.abp_pnorm{w}(p),abp_norm(abp_tnorm(p,:) == wtbl.abp_pnorm{w}(p)),'o','MarkerSize',10,'color',[0.08,0.00,0.5],'Linewidth',2)
% 
%             ylim([0 1.1])
%         end
%            %graph_beautify(gca,30,0.5,1)
% 
%     end
end
%Bin and average
%% PULSE AVERAGING
wtbl.bins_ln = placehold;
wtbl.tbin_ln = placehold;
wtbl.avg_bfi_ln = placehold;
wtbl.avg_abp_ln = placehold;
wtbl.std_bfi_ln = placehold;
wtbl.std_abp_ln = placehold;
wtbl.t_dicron_ln = placehold;
wtbl.loc_dicron_ln = placehold;
wtbl.t_next_ln = placehold;
wtbl.loc_next_ln = placehold;
wtbl.psf_ln = placehold;
wtbl.psp_ln = placehold;
wtbl.psf_loc_ln = placehold;
wtbl.psp_loc_ln = placehold;
wtbl.p_auc_ln = placehold;
wtbl.f_auc_ln = placehold;
wtbl.binsize_ln = placehold;
wtbl.f_auc_subedf = placehold;
wtbl.p_auc_subedp = placehold;

% Bins are centered on [0 1] but encompass BFI points slightly left shifted
% from 0
for w=1:height(wtbl)
    bins = -(bin_size/2):bin_size:(1+(bin_size/2));
    nbins = length(bins);

    bfi_tnorm = wtbl.bfi_tnorm{w};
    if isempty(bfi_tnorm)
        continue
    end
    abp_tnorm = wtbl.abp_tnorm{w};
    bfi_pstack = wtbl.bfi_pstack{w};
    abp_pstack = wtbl.abp_pstack{w};
    % bfi_clr = wtbl.clr{w}(:);
    % state = wtbl.substate{w};
    % sds_p = wtbl.sds(w);
    % name = wtbl.name{w};
    %wnum = wtbl.wnum(w);
    
    t_dicron = mean(wtbl.abp_dnorm{w});
    %t_next = mean(wtbl.abp_pnorm{w});
    t_next = median(wtbl.abp_pnorm{w});
    bin_dicron = zeros(length(bins)-1,1);
    bin_next = zeros(length(bins)-1,1);
    
    avg_bfi_ln = [];
    avg_abp_ln = [];
    std_bfi_ln = [];
    std_abp_ln = [];
    % if w==1
    %     xline(bins(1))
    % 
    % end
    % Binning ABP and BFI pulses
    for b = 2:nbins
        bfi_bin_ln = bfi_tnorm <= bins(b) & bfi_tnorm > bins(b-1);
        if t_dicron <= bins(b) & t_dicron > bins(b-1);
            bin_dicron(b-1) = 1;
        end
        if t_next <= bins(b) & t_next > bins(b-1);
            bin_next(b-1) = 1;
        end        
        avg_bfi_ln(b-1) = mean(bfi_pstack(bfi_bin_ln),'omitnan');
        std_bfi_ln(b-1) = std(bfi_pstack(bfi_bin_ln),'omitnan');
        abp_bin_ln = abp_tnorm <= bins(b) & abp_tnorm > bins(b-1);
        avg_abp_ln(b-1) = mean(abp_pstack(abp_bin_ln),'omitnan');
        std_abp_ln(b-1) = std(abp_pstack(abp_bin_ln),'omitnan');
       % if w==1
       % 
       %  xline(bins(b))
       %  plot(abp_tnorm',(abp_pstack./mean(mean(abp_pstack,'omitnan')))','red')
       %  plot(bfi_tnorm',(bfi_pstack./mean(mean(bfi_pstack,'omitnan')))','blue')
       % 
       % end
    end
    
%     if w==1
%         plot(bins,avg_abp_ln./(mean(avg_abp_ln)))
%     end

%% MOROPHOLOGICAL QUANTIFICATION
    % Find peak systolic flow and pressure of the resulting average pulse
    [fp,floc]=findpeaks(avg_bfi_ln);
    [pp,ploc]=findpeaks(avg_abp_ln);
    if isempty(floc)
        wtbl.psf_ln(w) = {NaN};
        wtbl.psf_loc_ln(w) = {NaN};
    else
        wtbl.psf_ln(w) = {fp(1)};
        wtbl.psf_loc_ln(w) = {floc(1)};
    end
        
    
    wtbl.psp_ln(w) = {pp(1)};
    
    wtbl.psp_loc_ln(w) = {ploc(1)};
    
    wtbl.binsize_ln(w) = {bin_size};
    wtbl.bins_ln(w) = {bins};
    wtbl.tbin_ln(w) = {bins(1:end-1) + bin_size/2};
    wtbl.t_dicron_ln(w) = {wtbl.tbin_ln{w}(logical(bin_dicron))};
    wtbl.loc_dicron_ln(w) = {find(bin_dicron)};
    wtbl.t_next_ln(w) = {wtbl.tbin_ln{w}(logical(bin_next))};
    wtbl.loc_next_ln(w) = {find(bin_next)};    
    wtbl.avg_bfi_ln(w) = {avg_bfi_ln};
    wtbl.std_bfi_ln(w) = {std_bfi_ln};
    wtbl.avg_abp_ln(w) = {avg_abp_ln};
    wtbl.std_abp_ln(w) = {std_abp_ln};
    wtbl.p_auc_ln(w) = {trapz(wtbl.tbin_ln{w},avg_abp_ln)};
    wtbl.f_auc_ln(w) = {trapz(wtbl.tbin_ln{w},avg_bfi_ln)};
    % area under the curve
    wtbl.f_auc_subedf(w) = {trapz(wtbl.tbin_ln{w},avg_bfi_ln-avg_bfi_ln(end))};
    wtbl.p_auc_subedp(w) = {trapz(wtbl.tbin_ln{w},avg_abp_ln-avg_abp_ln(end))};
end

gp = ~(cellfun(@isnan,wtbl.abp_avg_l));
wtbl.osf_ln = placehold;
wtbl.edf_ln = placehold;
wtbl.mf_ln = placehold;
wtbl.osp_ln = placehold;
wtbl.edp_ln = placehold;
wtbl.mp_ln = placehold;
wtbl.fph_ln = placehold;
wtbl.pph_ln = placehold;
wtbl.RI_ln = placehold;
wtbl.PI_ln = placehold;
wtbl.starts_ln = placehold;
wtbl.peaks_ln = placehold;

% Onset systolic flow, start of each pulse
wtbl.osf_ln(gp) = num2cell(cellfun(@(x) x(1), wtbl.avg_bfi_ln(gp)));
% End diastolic flow, end of each pulse
wtbl.edf_ln(gp) = num2cell(cellfun(@(x) x(end), wtbl.avg_bfi_ln(gp)));
% Mean of each pulse
wtbl.mf_ln(gp) = num2cell(cellfun(@(x,y) x./(bin_size*length(y)), wtbl.f_auc_ln(gp),wtbl.tbin_ln(gp)));

% Flow peak height (i.e. amplitude)
wtbl.fph_ln(gp) = num2cell(cell2mat(wtbl.psf_ln(gp)) - cell2mat(wtbl.edf_ln(gp)));
% Pulsatility and resisitivity indices
wtbl.PI_ln(gp) = num2cell(cell2mat(wtbl.fph_ln(gp))./cell2mat(wtbl.mf_ln(gp)));
wtbl.RI_ln(gp) = num2cell(cell2mat(wtbl.fph_ln(gp))./cell2mat(wtbl.psf_ln(gp)));

% Same morphological params as above but for ABP
wtbl.osp_ln(gp) = num2cell(cellfun(@(x) x(1), wtbl.avg_abp_ln(gp)));
wtbl.edp_ln(gp) = num2cell(cellfun(@(x) x(end), wtbl.avg_abp_ln(gp)));
wtbl.mp_ln(gp) = num2cell(cellfun(@(x,y) x./(bin_size*length(y)), wtbl.p_auc_ln(gp),wtbl.tbin_ln(gp)));
wtbl.pph_ln(gp) = num2cell((cell2mat(wtbl.psp_ln(gp)) - cell2mat(wtbl.edp_ln(gp))));

cnames = categorical(wtbl.name);
names = categories(cnames);
binsdiff = diff(wtbl.bins_ln{1});
bin = binsdiff(1);

wtbl.starts_ln(gp) = num2cell(cellfun(@(x) mean(x(:,1)),wtbl.bfi_tnorm(gp)));
wtbl.peaks_ln(gp) = num2cell(cell2mat(wtbl.psf_loc_ln(gp)).*bin);

%% Dat Tbl - Data sorted by subject rather than one row for each window, to make things more manageable
substates = unique(wtbl.substate)
dat_tbl = table();
dat_tbl.name = unique(wtbl.name); % assumes we are working with a single subject at a time
for ss=1:length(substates)
    sub_tbl = wtbl(categorical(wtbl.substate) == substates{ss},:);

    %%%%%%%% Window-averaged pulse QC %%%%%%%%%%
    sub_tbl_qc = pulse_tests3_1(sub_tbl);
    %sub_tbl_qc = pulse_tests4(sub_tbl);
    %sub_tbl_qc = pulse_tests3(sub_tbl);
    dat_tbl(1,substates{ss}) = {sub_tbl_qc};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end