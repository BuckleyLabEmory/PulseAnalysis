function state_pulse_plot(dat_tbl,state,save_path)

%%
dat_tbl.(state)
%state_name = dat_tbl.(state){1}.state{1};
all_dat = vertcat(dat_tbl.(state){:});
% IDs no identifiers
ids = string(cellfun(@(x) strsplit(x,' '),all_dat.name,'UniformOutput',false));
subs = categorical(ids);
sub_i = categories(subs);

    n_subs = length(sub_i);

    
    for n=1:length(sub_i)
        f=figure
        set(f,'defaultAxesColorOrder',[[0 0 1]; [0.6 0 0]])
        hold on
        
        yyaxis left
        tbl = all_dat(subs == sub_i(n),:);
        pt = tbl.tbin_ln{1};
        fpulses = vertcat(tbl.avg_bfi_ln{:});
        fpm = median(fpulses,1);
        fpstd = std(fpulses);

        pulse_fillplot(pt,fpm,fpstd,'blue');
        graph_beautify(gca,15,2,2)
        hold on 
        
        yyaxis right
        ppulses = vertcat(tbl.avg_abp_ln{:});
        ppm = median(ppulses,1);
        ppstd = std(ppulses);
        clr = [0.6 0 0];

        pulse_fillplot(pt,ppm,ppstd,clr);
        title([{['Subject: ' sub_i{n}]}...
            {['state: ' state]}...
            {'average pulses abp and bfi'}],'Interpreter','none')
        graph_beautify(gca,15,2,2)
    end
    %sgtitle([string(state) "Subject average pulses abp and bfi"])
    id = [state ' abp and bfi']
    if ~isempty(save_path)
        saveas(gca,[save_path filesep id ' pulses.png'],'png')
        savefig(gcf,[save_path filesep id ' pulses.fig'])
    end
end