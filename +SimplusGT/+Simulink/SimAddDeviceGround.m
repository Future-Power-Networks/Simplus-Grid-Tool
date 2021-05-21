% This function adds apparatus ground

% Author(s): Yitong

function SimAddDeviceGround(Name_Model,Size_D_GND,Shift_D_GND,FullName_Device,Name_Device,DeviceType)

% Organize data
N_Device = length(DeviceType);

% Add blocks
for i = 1:N_Device

  	% If the apparatus is an active apparatus and an ac apparatus
    if DeviceType{i} <= 90
        
        % Add apparatus ground
        Name_D_GND{i} = ['D-GND' num2str(i)];
        FullName_D_GND{i} = [Name_Model '/' Name_D_GND{i}];
        add_block('powerlib/Elements/Ground',FullName_D_GND{i});
        PortPos_Device{i} = get_param(FullName_Device{i},'PortConnectivity');
        % Position of apparatus ground is set by referring to the position of
        % corresponding apparatus
        Pos_D_GND{i} = PortPos_Device{i}(5).Position;
        Pos_D_GND{i} = Pos_D_GND{i} + Shift_D_GND;
        set_param(FullName_D_GND{i},'position',[Pos_D_GND{i},Pos_D_GND{i}+Size_D_GND]);
        
        % Connect apparatus to apparatus ground
        add_line(Name_Model,[Name_Device{i} '/LConn4'],[Name_D_GND{i} '/LConn1'],...
                'autorouting','smart');
        
    end

end

end