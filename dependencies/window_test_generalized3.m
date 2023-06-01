function [tests,result,window_BFI_update,window_ABP_update,window_keep,discard_reason] = window_test_generalized3(window_BFI,window_ABP,window_time,numnans)
tests = [];
result = [];
test_count = 0;
window_keep = 1;
discard_reason = 'none';
window_BFI_update = window_BFI;
window_ABP_update = window_ABP;



%% NaNs in both forward and backward shift next idx
% test_result = strfind(tests,'Both first and last NaN');
%     if result(test_result) == 1
%         test_count = test_count + 1;
%         tests{test_count} = 'Backward or forward shift by 1 both NaNs';
%         fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
% 
%         lower_forwardupdate = pdo.icp_lower(l)-1;
%         upper_forwardupdate = pdo.icp_upper(l)-1;
%         lower_backwardupdate = pdo.icp_lower(l)+1;
%         upper_backwardupdate = pdo.icp_lower(l)+1;
% 
% 
%         if all(isnan([lower_forwardupdate;upper_forwardupdate;lower_backwardupdate;upper_backwardupdate]))
%            result(test_count) = 1;
%            fprintf('Result: True\n')
%            window_keep = 0;
%            return
%         else
%             result(test_count) = 0;
%             fprintf('Result: False\n')
%         end
%         
%         if isnan(window_BFI(1))
%             forward_shift = 1;
%             backward_shift = 0;
%         elseif 
%     end
% 
% %% Forward shift possible    
% test_result = strfind(tests,'Both first and last NaN');
%     if result(test_result) == 1    
%         % NaN at first position
%         test_count = test_count + 1;
%         tests{test_count} = 'Shift forward';
%         fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
% 
%             if isnan(window_BFI(1))
%                pdo_icp_lower_update = pdo.icp_lower(l)-1;
%                pdo_icp_upper_update = pdo.icp_upper(l)-1;
%                 
%                pdo.icp_lower(l) = pdo_icp_lower_update;
%                pdo.icp_upper(l) = pdo_icp_upper_update;
%                result(1) = 0;
%                window_BFI_update = DCS_resample(pdo.icp_lower(l):pdo.icp_upper(l));
%                window_ABP_update = DCS
%             end
%     end
% %% Backward shift possible
% 
% %% NaN at last position
% 
%     if isnan(window_BFI(end))
%         if result(1) == 0
%             
%         end
%         pdo_icp_lower_update = pdo.icp_lower(l)+1;
%         pdo.icp_lower(1) = pdo_icp_lower_update;
%         pdo_icp_upper_update = pdo.icp_upper(l)+1;
%         pdo.icp_upper(l) = pdo_icp_upper_update;
%         result(2) = 0;
%     end
%% Multiple NaN Test 
    test_count = test_count + 1;
    tests{test_count} = 'All NaNs';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    if all(isnan(window_BFI))
        result(test_count) = 1;
        fprintf('Result: True\n')
        window_keep = 0;
        discard_reason = 'All NaNs';
        return
    end
    
%% Multiple NaN Test 
    test_count = test_count + 1;
    tests{test_count} = 'Sequential NaNs';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
%     % Test ABP
    [~,~,aIB] = RunLength(window_ABP);
    t = diff(aIB)>numnans;
    if any(t)
        result(test_count) = 1;
        fprintf('Result: True\n')
        window_keep = 0;
        discard_reason = 'NaNs Threshold - ABP';
        return
    else
    result(test_count) = 0;
    fprintf('Result: False\n')
    end
    
    test_count = test_count + 1;
    tests{test_count} = 'Sequential NaNs';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    % Test BFI
    [~,~,bIB] = RunLength(window_BFI);
    t = diff(bIB)>numnans;
    if any(t)
        result(test_count) = 1;
        fprintf('Result: True\n')
        window_keep = 0;
        discard_reason = 'NaNs Threshold - BFI';
        return
    else
        result(test_count) = 0;
        fprintf('Result: False\n')
    end
%     


%% NaN at both first or last position
    test_count = test_count + 1;
    tests{test_count} = 'Both first and last NaN';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    % Test BFI
    if isnan(window_BFI(1)) | isnan(window_BFI(end))
       result(test_count) = 1;
       fprintf('Result: True\n')
       window_keep = 0;
       discard_reason = 'NaN in first or last idx';
       %return
    else
        result(test_count) = 0;
        window_BFI_update = window_BFI;
        fprintf('Result: False\n')
    end
    
    % Test ABP
    if isnan(window_ABP(1)) | isnan(window_ABP(end))
       result(test_count) = 1;
       fprintf('Result: True\n')
       window_keep = 0;
       discard_reason = 'NaN in first or last idx';
       %return
    else
        result(test_count) = 0;
        window_ABP_update = window_ABP;
        fprintf('Result: False\n')
    end    
    
%% Any NaNs 
    
    test_count = test_count + 1;
    tests{test_count} = 'Any NaNs';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])   
    
    % Test BFI
    if any(isnan(window_BFI))
        result(test_count) = 1;
        fprintf('Result: True. Interpolating...\n')
        window_BFI_update = interp_nans(window_BFI,window_time,0);
    else
        result(test_count) = 0;
        fprintf('Result: False \n')
    end

    % Test APB
    if any(isnan(window_ABP))
        result(test_count) = 1;
        fprintf('Result: True - ABP. Interpolating...\n')
        window_ABP_update = interp_nans(window_ABP,window_time,0);
    else
        result(test_count) = 0;
        fprintf('Result: False \n')
    end
%% CNAP calibration
    test_count = test_count + 1;
    tests{test_count} = 'CNAP calibration';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])

    if length(find(window_ABP < 5)) > 10
        result(test_count) = 1;
        fprintf('Result: True - CNAP calibration...\n')
        window_ABP_update = interp_nans(window_ABP,window_time,0);
    else
        result(test_count) = 0;
        fprintf('Result: False \n')
    end
    

%% Noise characteristics test

% if ~any(result)
%     window_BFI_update = window_BFI;
% 
% end
end