% This function adds scopes for apparatuses.

% Author(s): Yitong Li

function SimAddDeviceScope(Name_Model,Size_D_Scope,Shift_D_Scope,Size_DS_Bus,Shift_DS_Bus,Pos_Device,Name_Device,DeviceType)

% Organize data
N_Device = length(DeviceType);      % Number of apparatuses

% Add block
for i = 1:N_Device
    
    % If the device is not an floating bus
    if DeviceType{i}~=0100 && DeviceType{i}~=1100

        % Add device scope bus
        Name_DS_Bus{i} = ['DS-Bus' num2str(i)];
        FullName_DS_Bus{i} = [Name_Model '/' Name_DS_Bus{i}];
        add_block('simulink/Signal Routing/Bus Selector',FullName_DS_Bus{i});
        set_param(gcb,'Orientation','left');
        Pos_DS_Bus{i} = Pos_Device{i} + Shift_DS_Bus;
        set_param(gcb,'position',[Pos_DS_Bus{i},Pos_DS_Bus{i}+Size_DS_Bus]);
        if DeviceType{i}<=90
            Output_DS_Bus = ['v_dq,i_dq,v_abc,i_abc,w,theta,pq'];
            Length_D_Measurement = 7;
        elseif 1000<=DeviceType{i} && DeviceType{i}<=1090
            Output_DS_Bus = ['v,i,p'];
            Length_D_Measurement = 3;
        elseif 2000<=DeviceType{i} && DeviceType{i}<=2090
            Output_DS_Bus = ['v_dq,i_dq,v_abc,i_abc,w,theta,pq_ac,v,i,p_dc'];
            Length_D_Measurement = 10;
        end
        
        set_param(gcb,'OutputSignals',Output_DS_Bus);
        Port_DSBus{i} = get_param(gcb,'PortHandles');

        % Conect scope bus to device
        add_line(Name_Model, {[Name_Device{i} '/1']}, {[Name_DS_Bus{i} '/1']});

        % Add device scope
        Name_D_Scope{i} = ['D-Scope' num2str(i)];
        FullName_D_Scope{i} = [Name_Model '/' Name_D_Scope{i}];
        add_block('simulink/Sinks/Scope',FullName_D_Scope{i});
        set_param(gcb,'NumInputPorts',num2str(Length_D_Measurement));
        set_param(gcb,'LayoutDimensionsString',['[' num2str(Length_D_Measurement) ' 1]']);
        Port_D_Scope{i} = get_param(gcb,'PortHandles');
        set_param(gcb,'Orientation','left');
        Pos_D_Scope{i} = Pos_Device{i} + Shift_D_Scope;
        set_param(gcb,'position',[Pos_D_Scope{i},Pos_D_Scope{i}+Size_D_Scope]);

        % Connect scope to scope bus
        for k = 1:Length_D_Measurement
            add_line(Name_Model,{Port_DSBus{i}.Outport(k)},{Port_D_Scope{i}.Inport(k)});
        end
        
    end
end

end