function tbl_qc = pulse_tests3_1(tbl)
%% This version differs from pulse_tests2 in that it just focuses on sds2, because for the 2022 CVR pulsatility manuscript that is the most important thing. Rejected ABP pulses lead to overall rejection still.
    % for wi = 1:height(tbl)
    %     counts(wi) = sum(tbl.wnum==tbl.wnum(wi)); % number of times each unique value is repeated  
    % end
    % r_pulses = counts > 1;
    % tbl1 = tbl(r_pulses,:);

%%
    tbl1 = tbl;
    winds = unique(tbl1.wnum);   
    tbl_qc = tbl1;
    tbl_qc.wp_accept = cell(height(tbl_qc),1);
%%    

    for w=1:length(winds)
    % pulses are tested by window - if one sds or abp pulse fails they are all be rejected to keep the dataset paired    
        
        % Only keeping pulses for which we have both SDS1 and 2

        
        wselect = tbl1.wnum==winds(w);
        tbl_rw = tbl1(wselect,:);
        
        tpulse = tbl_rw.tbin_ln{1}; % the same across pulses and windows
        
        pulses = {};
        try
%         pulses{1} = tbl_rw.avg_bfi_ln{1}; % sds 1 pulse from window
%         pulses{2} = tbl_rw.avg_bfi_ln{2}; % sds2 pulse from window
%         %pulses{3} = tbl_rw.avg_abp_ln{1}; % abp pulse from window
        pulses{1} = tbl_rw.avg_bfi_ln{1}; % sds2 pulse from window
        pulses{2} = tbl_rw.avg_abp_ln{1}; % abp pulse from window
        %pulses{3} = tbl_rw.avg_abp_ln{1}; % abp pulse from window
        catch
            fprintf("Error")
        end
        pulse_peak = {};
        pulse_peak{1} = tbl_rw.psf_ln{1};
        pulse_peak{2} = tbl_rw.psp_ln{1};
 
        pulse_height = {};
        pulse_height{1} = tbl_rw.fph_ln{1};
        pulse_height{2} = tbl_rw.pph_ln{1};

        pulse_end = {};
        pulse_end{1} = tbl_rw.edf_ln{1};
        pulse_end{2} = tbl_rw.edp_ln{1};
        
        if tbl_rw.avg_abp_ln{1} ~= tbl_rw.avg_abp_ln{1} % this shouldn't happen
            pdiff = tbl_rw.avg_abp_ln{1} - tbl_rw.avg_abp_ln{2};
            warning(['Discovered abp pulse discrepancy of ' num2str(mean(pdiff)) ' mmHg avg'])
        end
        
        accept = ones(1,length(pulses));
        
        for p=1:length(pulses)
            
            pulse = pulses{p};
            ph = pulse_height{p};
            pp = pulse_peak{p};
            pe = pulse_end{p};
            t_pp = tpulse(pulse == pp);
            t_pe = tpulse(pulse == pe);
            if find(pulse == pe) ~= length(pulse)
                error('edf not last timepoint in pulse')
            end
            
            
             % test 1: are there spikes larger than the systolic peak later in the
             % pulse?
            end_pulse = pulse(tpulse>0.6);
            if any(end_pulse > max(pulse))
                accept(p) = 0;
            else 
            end
 
            % test 2: is the pulse height negative?
            if ph < 0 
                accept(p) = 0;
            else
            end

            % test 3: If the slope of systolic period of pulse is negative
            pulse_start = pulse(tpulse<t_pp);
            if length(pulse_start) < 3
               upstrk = pulse_start;
            else 
                upstrk = pulse_start(1:3);
            end
            if mean(diff(upstrk)) < 0
%             if mean(pulse_start - pulse(1)) < 0
                accept(p) = 0;
            end
            
            %error('Look into this and make sure it is doing what I think it is doing')
            if diff(pulse_start - pulse(1)) < 0
                accept(p) = 0;
            end

            % test 4: If end diastolic flow varies by more than 200% from
            % onset systolic flow
            if pulse(end)/pulse(1) > 2
                accept(p) = 0;
            end

            % test 6: is the edf greater than the mean of the pulse
            if pe > mean(pulse)
                accept(p) = 0;
            end            
          
            % test 7: Is the ABP pulse is very low from ABP cal
            if p==2 & all(pulse<10)
                accept(p) = 0;
            end
            
        end
        
        if any(accept == 0)
            %tbl_qc.wp_accept(wselect) = {0,0};   
            tbl_qc.passed_pulse_tests(wselect) = 0;
        elseif all(accept == 1)
            %tbl_qc.wp_accept(wselect) = {1,1};
            tbl_qc.passed_pulse_tests(wselect) = 1;
        end 
    end
end