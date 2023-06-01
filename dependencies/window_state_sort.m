function wtbl_n = window_state_sort(ftbl,window_start,window_end)
% This code associates each window with a state and substate based on the
% given marks table. It has not been tested thoroughly for inputs formatted
% differently than the examplies i.e. with multiple states or with more than two substates

    wstart = ftbl.t_bfi{:}(window_start(window_end < length(ftbl.t_bfi{:})));
    wend = ftbl.t_bfi{:}(window_end(window_end < length(ftbl.t_bfi{:})));
    wbounds = cell2table([num2cell(wstart),num2cell(wend)]);
    wbounds.Properties.VariableNames = {'wstart','wend'};
    %winfo = repmat(ftbl(1,{'hv','sess','name'}),length(wend),1);
    winfo = repmat(ftbl(1,{'name','state'}),length(wend),1);
    wtbl_n = [winfo, wbounds];
    % Initialize state as "none". Will stay like this if no marks given.
    startstate = cellstr(repmat("none",length(wend),1)); 

    wtbl_n.state = startstate;
    wtbl_n.substate = startstate;

    if any(string(ftbl.Properties.VariableNames)=="run_tbl")
        run_tbl = ftbl.run_tbl{1};
        marks_tbl = ftbl.marks{1};
        varnames = run_tbl.Properties.VariableNames';
        statenames = cellfun(@any, cellfun(@(x) contains(x,'State','IgnoreCase',true),varnames,'UniformOutput',false));
        instate = varnames(statenames);
        issubstate = cellfun(@any, cellfun(@(x) contains(x,'substate','IgnoreCase',true),instate,'UniformOutput',false));
        isstate = ~issubstate;
        states = instate(isstate);
        state_loc = cellfun(@(x) strfind(x,'State'),states,'UniformOutput',false);
        state_names = cellfun(@(x,y) x(1:y-1),states,state_loc,'UniformOutput',false);
        state_names = unique(state_names);
        state_starts = cellfun(@any,cellfun(@(x) contains(x,'StateStart'),instate,'UniformOutput',false));
        state_ends = cellfun(@any,cellfun(@(x) contains(x,'StateEnd'),instate,'UniformOutput',false));
        substate_starts = cellfun(@any,cellfun(@(x) contains(x,'SubstateStart'),instate,'UniformOutput',false));
        substate_ends = cellfun(@any,cellfun(@(x) contains(x,'SubstateEnd'),instate,'UniformOutput',false));
        instate(isstate & state_starts)
        instate(isstate & state_ends)
    
        sname = {};
        sn = 1;
        for s=1:length(state_names)
            state_names{s};
            all_in = cellfun(@(x) contains(x,state_names{s}),instate);
            start_mark = run_tbl.(instate{all_in & state_starts});
            end_mark = run_tbl.(instate{all_in & state_ends});
            start_time = marks_tbl(marks_tbl(:,2) == start_mark,3);
            end_time = marks_tbl(marks_tbl(:,2) == end_mark,3);
            sname{sn,1} = state_names{s};
            sname{sn,2} = [start_time,end_time];
            state_substates = instate(all_in & issubstate & substate_starts);
            substate_loc = cellfun(@(x) strfind(x,'Substate'),state_substates,'UniformOutput',false);
            state_substate_names = cellfun(@(x,y) x(1:y-1),state_substates,substate_loc,'UniformOutput',false);
            if ~isempty(state_substates)
            else
                sn = sn+1;
            end
    
            for ss = 1:length(state_substate_names)
                all_in_substate = cellfun(@(x) contains(x,state_substate_names{ss}),instate);
                substate_start_mark = run_tbl.(instate{all_in_substate & substate_starts});
                substate_start_time = marks_tbl(marks_tbl(:,2) == substate_start_mark,3);
                substate_end_mark = run_tbl.(instate{all_in_substate & substate_ends});
                substate_end_time = marks_tbl(marks_tbl(:,2) == substate_end_mark,3);
                sname{sn,1} = state_names{s};
                sname{sn,2} = [start_time,end_time];
                sname{sn,3} = state_substate_names{ss};
                sname{sn,4} = [substate_start_time, substate_end_time];
                sn = sn + 1;    
            end
    
        end
    
        sts = unique(sname(:,1));
        for st = 1:size(sts,1) % for each state
            cstate = sts{st}
            cstate_array = find(cell2mat(cellfun(@(x) contains(x,cstate),sname(:,1),'UniformOutput',false)));
            sbounds = sname{cstate_array(1),2};
            if isempty(sbounds) % if this session doesn't have this state
                continue
            end
            is_instate = sbounds(1) <= wtbl_n.wstart & wtbl_n.wend <= sbounds(2);
            wtbl_n(is_instate,"state") = {cstate};
            if size(sname,2) < 3 % if there were no substates these fields won't be here and the below loop is irrelevant
                continue
            else
                for sb = 1:length(cstate_array)
                        csubstate = sname{cstate_array(sb),3} 
                    if isempty(csubstate)
                        continue
                    end
                    sb_bounds = sname{cstate_array(sb),4};
                    is_in_substate = sb_bounds(1) <= wtbl_n.wstart & wtbl_n.wend <= sb_bounds(2);
                    wtbl_n(is_in_substate,"substate") = {csubstate};
                end
            end
        end
    end
end