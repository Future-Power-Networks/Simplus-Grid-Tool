function DeviceMaskInitialization(EnableInsideModification,DiscreDampingFlag)
    if EnableInsideModification == 1
        if DiscreDampingFlag == 0
            set_param([gcb '/R_VD'],'Commented','on');
        else
            set_param([gcb '/R_VD'],'Commented','off');
        end
    %     if DiscreDampingFlag == 0
    %         delete_block([gcb '/R_VD']);
    %     end
    end
end