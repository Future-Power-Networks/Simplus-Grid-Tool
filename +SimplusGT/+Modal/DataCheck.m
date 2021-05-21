function DataCheck(AxisSel, ApparatusSelL12, ModeSelAll, ApparatusSelL3All,...
    StateSel_DSS, ModeSel_DSS,BodeEnable,Layer12Enable,Layer3Enable,StatePFEnable)
if StatePFEnable == 1
    if ModeSel_DSS == 0
        error('Please select at least a mode for State-PF analysis');
    elseif StateSel_DSS == 0
        error('Please select at least a state for State-PF analysis');
    end
end 

if BodeEnable == 1
    if AxisSel == 0
        error('Please select at least an axis for bode plot');
    elseif ApparatusSelL12 == 0
        error('Please select at least a apparatus for bode plot');
    end
end

if Layer12Enable == 1
    if ApparatusSelL12 == 0
        error('Please select at least a an apparatus for Impedance-PF Layer 1&2 analysis');
    elseif ModeSelAll == 0
        error('Please select at least a mode for for Impedance-PF Layer 1&2 analysis');
    end
end

if Layer3Enable == 1
    if ApparatusSelL3All == 0
        error('Please select at least a an apparatus for Impedance-PF Layer 3 analysis');
    elseif ModeSelAll == 0
        error('Please select at least a mode for for Impedance-PF Layer 3 analysis');
    end
end
    
    
end