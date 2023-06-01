% generate data table
clear
close all
path = '/Volumes/labs/buckley-lab/Projects/Waveform_analysis/0_Papers/2023_PulsatilityPaper/40_PULSEANALYSIS_final/'
load([path 'sub_stats_bfi.mat']);
load([path 'sub_stats_abp.mat'])
load([path 'fd_td_crcp.mat'])
load([path 'dat_tbl_qc.mat'])
load([path 'crcp_tbl_qc.mat'])
%%
cohort = readtable('/Volumes/labs/buckley-lab/Projects/Waveform_analysis/0_Papers/2023_PulsatilityPaper/10_DATA_LOADED_final/FinalCohort_DCS.xlsx')
%%
savepath = '/Volumes/labs/buckley-lab/Projects/Waveform_analysis/0_Papers/2023_PulsatilityPaper/40_PULSEANALYSIS_final/'
%%

abp_params = {
            'onset pressure',               'osp_ln';...
            'peak systolic pressure',       'psp_ln';...
            'end diastolic pressure',       'edp_ln';...            
            'pressure peak height',         'pph_ln';...
            'pressure under the curve',     'p_auc_subedp';...
            'mean pressure',                'mp_ln';...  
};

bfi_params = {
            'onset flow',                   'osf_ln';...
            'peak systolic flow',           'psf_ln';...
            'end diastolic flow',           'edf_ln';...
            'flow peak height',             'fph_ln';... 
            'flow area under the curve',    'f_auc_subedf';...
            'mean flow',                    'mf_ln';...
            'flow resistive index',         'RI_ln';...
            'flow pulsatility index',       'PI_ln';...

};

crcp_params = {
            'TD CrCP - runoff',             'td_crcp';...
            'TD CrCP - whole pulse',        'td_crcp_pulse';...
            'FD CrCP',                      'fd_crcp';...
      };

all_params = [bfi_params ; abp_params ; crcp_params];
n_params = size(all_params,1);

%%
abp = sub_stats_abp(sub_stats_abp.sds==2 & ismember(string(sub_stats_abp.param),string(abp_params(:,2)))...
    ,{'name','state_name','param','n_used','mean','sdev'});
bfi = sub_stats_bfi(sub_stats_bfi.sds==2 & ismember(string(sub_stats_bfi.param),string(bfi_params(:,2)))...
    ,{'name','state_name','param','n_used','mean','sdev'});
tmp = vertcat(vertcat(crcp_tbl_qc.hc_bl{:}),vertcat(crcp_tbl_qc.hc_end{:}));
crcp_all = tmp(tmp.sds==2,[{'hv','substate'} crcp_params(:,2)'])
crcp_all.hv = arrayfun(@(x) strjoin(strsplit(x,"_")," "), string(crcp_all.hv));
crcp_proc.hv = arrayfun(@(x) strjoin(strsplit(x,"_")," "), string(crcp_proc.hv));
%%
subs = unique(crcp_proc.hv);
states = unique(crcp_all.substate);
state_names = {"Baseline","Hypercapnia"};
crcp = table();
for p=1:size(crcp_params,2) % for every parameter, make a column vector and add to the table
    crcp_p = crcp_proc(:,[{'hv','substate'} crcp_params{p,2}])
    param_col = repmat(string(crcp_params{p,2}),height(crcp_p),1)
    crcp_p.param = param_col;
    crcp_p.state_name = cell(height(crcp_p),1);
    crcp_p.n_used = zeros(height(crcp_p),1);
    crcp_p.mean = zeros(height(crcp_p),1);
    crcp_p.sdev = cell(height(crcp_p),1);

    for s=1:length(subs) % then populate mean and standard deviation for each subject. Mean needs to come from CrCP proc because of twa()
        for ss=1:length(states) 
            row_i = categorical(crcp_p.substate) == states{ss} & categorical(crcp_p.hv) == subs{s};
            sdev_rows = categorical(crcp_all.substate) == states{ss} & categorical(crcp_all.hv) == subs{s};
            sdev_dat = crcp_all{sdev_rows,crcp_params{p,2}};
            crcp_p.n_used(row_i) = length(sdev_dat(~sdev_dat == 0));
            crcp_p.mean(row_i) = crcp_p{row_i,crcp_params{p,2}};
            crcp_p.sdev{row_i} = std(sdev_dat(~sdev_dat == 0));
            crcp_p.state_name{row_i} = state_names{ss};
        end
    end
    crcp = [crcp;crcp_p(:,{'hv','state_name','param','n_used','mean','sdev'})];
end
% Just rename "hv" to be "name"
crcp.Properties.VariableNames{1} = 'name';
%% all_dat
all_dat = [abp;bfi;crcp];
%% Add age and sex data
all_dat.age = zeros(height(all_dat),1);
all_dat.sex = cell(height(all_dat),1);
for s=1:length(subs)
    fill_rows = categorical(all_dat.name) == subs{s};
    sub_info = cohort(categorical(cohort.ID)==subs{s},:);
    all_dat.age(fill_rows) = repmat(sub_info.age,length(find(fill_rows)),1);
    all_dat.sex(fill_rows) = repmat(sub_info.sex,length(find(fill_rows)),1);
end
%% Rename parameters to their full names
all_dat2 = all_dat;
for p=1:n_params
    pfind = all_params{p,2};
    prows = categorical(all_dat.param) == pfind;
    all_dat2.param(prows) = repmat(all_params{p,1},length(find(prows)),1);
end

all_dat2 = sortrows(all_dat2,'name')

%% Remove initial from ID
all_dat3 = all_dat2;
all_dat3.name = cellfun(@(x) x(length(x)-3:length(x)),all_dat2.name,'UniformOutput',false)
clear all_dat2
%%
all_dat4 = all_dat3;
all_dat4 = movevars(all_dat4,'sex','before','state_name')
all_dat4 = movevars(all_dat4,'age','before','state_name')
%% Write table
writetable(all_dat4,[savepath 'healthy-adult-pulsemorpho.csv'])
%% Generate CrCP Tables - had to regenerate to get crcp_tbl_qc
% states = {'hc_bl','hc_end'}
% state_names = {'resting','hypercapnia'};
% method = 'median';
% tshift = 0;
% 
% [crcp_tbl_qc,crcp_tbl_raw,qc_result] = td_fd_crcp_CNAP_v2(dat_tbl_qc,states,method,tshift); % Calculates FD and TD CrCP (peak method and dicrotic notch method)
% 
% % for p=1:size(crcp_params,2)
% %     crcp_p = crcp_proc(:,[{'hv','substate'} crcp_params{p,2}])
% % end

