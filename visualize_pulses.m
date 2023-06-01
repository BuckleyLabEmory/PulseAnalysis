states = string(dat_tbl.Properties.VariableNames(2:end));
dat_tbl_qc = dat_tbl(:,'name');
for s=1:length(states)
    dat_tbl_qc.(states{s}) = cell(height(dat_tbl),1);
    for subs=1:height(dat_tbl)
        % %%%%% Use quality control filters
        qc = dat_tbl{subs,states{s}}{1}.passed_pulse_tests == 1;
        dat_tbl_qc(subs,states{s}) = {dat_tbl{subs,states{s}}{1}(qc,:)};
    end
    state_pulse_plot(dat_tbl_qc,states{s},[])
end

