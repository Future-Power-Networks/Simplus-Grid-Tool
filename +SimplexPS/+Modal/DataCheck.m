function DataCheck(AxisSel, DeviceSelL12, ModeSelAll, DeviceSelL3All,...
    StateSel_DSS, ModeSel_DSS,BodeEnable,Layer12Enable,Layer3Enable,StatePFEnable)
if StatePFEnable == 1
    if ModeSel_DSS == 0
        error('Please select at least a mode for State-PF anallysis');
    elseif StateSel_DSS == 0
        error('Please select at least a state for State-PF anallysis');
    end
end 

if BodeEnable == 1
    if AxisSel == 0
        error('Please select at least an axis for bode plot');
    elseif DeviceSelL12 == 0
        error('Please select at least a device for bode plot');
    end
end

if Layer12Enable == 1
    if DeviceSelL12 == 0
        error('Please select at least a adevice for Impedance-PF Layer 1&2 analysis');
    elseif ModeSelAll == 0
        error('Please select at least a mode for for Impedance-PF Layer 1&2 analysis');
    end
end

if Layer3Enable == 1
    if DeviceSelL3All == 0
        error('Please select at least a adevice for Impedance-PF Layer 3 analysis');
    elseif ModeSelAll == 0
        error('Please select at least a mode for for Impedance-PF Layer 3 analysis');
    end
end
    
    
end