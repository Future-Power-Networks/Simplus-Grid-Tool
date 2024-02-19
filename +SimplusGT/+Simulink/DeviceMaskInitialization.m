function ApparatusMaskInitialization(EnableInsideModification,DiscreDampingFlag)
    % Check if Rv needs to be commented out
    % Notes: A same block twice not be commented twice. But every time when
    % runing the code or when clicking the "apply" button of the mask, the
    % initialization code will be run. Therefore, an additional flag
    % "EnableInsideModification" is added to avoid this problem.
    if EnableInsideModification == 1
        if DiscreDampingFlag == 0
            set_param([gcb '/R_VD'],'Commented','on');
        else
            set_param([gcb '/R_VD'],'Commented','off');
        end
    end
    
%     % Set the value of Rv
%     Rv = 1;
    
end