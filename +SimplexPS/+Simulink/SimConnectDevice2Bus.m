% This function connects device to bus

% Author(s): Yitong Li

function SimConnectDevice2Bus(Name_Model,Name_Bus,Name_Device,DeviceType)

% Organize data
N_Device = length(DeviceType);

% Add block
for i = 1:N_Device
    
    % If the device is NOT a "floating bus"
    if floor(DeviceType{i}/10) <= 9
        add_line(Name_Model,...
            {[Name_Device{i} '/Lconn1'],[Name_Device{i} '/Lconn2'],[Name_Device{i} '/Lconn3']},...
            {[Name_Bus{i} '/Lconn1'],[Name_Bus{i} '/Lconn2'],[Name_Bus{i} '/Lconn3']},...
            'autorouting','smart');
    end
    
end

end