% This function adds devices into simulink model.

% Author(s): Yitong Li

function [FullName_Device,Name_Device,Pos_Device] = SimAddDevice(Name_Model,Name_LibFile,Size_Device,Shift_Device,Pos_Bus,DeviceType,ListAdvance,PowerFlow)

% Organize data
DiscreMethod = ListAdvance(1);
LinearizationTimes = ListAdvance(2);
DiscreDampingFlag = ListAdvance(3);
DirectFeedthrough = ListAdvance(4);

N_Device = length(DeviceType);

% Add device
for i = 1:N_Device
    if floor(DeviceType{i}/10) <= 9
        
        switch floor(DeviceType{i}/10)
            case 0
                Name_Device{i} = ['SM' num2str(i)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Synchronous Machine (dq-Frame System Object)'],FullName_Device{i});
            case 1
                Name_Device{i} = ['VSI-PLL' num2str(i)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Grid-Following Voltage-Source Inverter (dq-Frame System Object)'],FullName_Device{i});
            case 2
                Name_Device{i} = ['VSI-Droop' num2str(i)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Grid-Forming Voltage-Source Inverter (dq-Frame System Object)'],FullName_Device{i});
            case 9
            	Name_Device{i} = ['Inf-Bus' num2str(i)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Infinite Bus'],FullName_Device{i});
            otherwise
                error(['Error']);
        end
        
        % Set position
       	% The position of device is set by referring to the position of correpsonding bus
        Pos_Device{i} = Pos_Bus{i} + Shift_Device;
        set_param(FullName_Device{i},'position',[Pos_Device{i},Pos_Device{i}+Size_Device]);
        set_param(FullName_Device{i},'Orientation','left');
        
        % Set common variables
      	set_param(gcb,'Sbase','Sbase');
        set_param(gcb,'Vbase','Vbase');
        set_param(gcb,'Wbase','Wbase');
        set_param(gcb,'Ts','Ts');
        
        % If the device is an "active device"
        if floor(DeviceType{i}/10) <= 5
            
            % Set system object parameters
            set_param(gcb,'DeviceType',['DeviceType{' num2str(i) '}']);
            set_param(gcb,'DevicePara',['DevicePara{' num2str(i) '}']);
            set_param(gcb,'PowerFlow',['PowerFlow{' num2str(i) '}']);
            set_param(gcb,'x0',['x_e{' num2str(i) '}']);
            set_param(gcb,'OtherInputs',['OtherInputs{' num2str(i) '}']);

            % Set discretization method
            switch DiscreMethod
                case 1
                    DeviceDiscreMethod = 'Forward Euler';
                case 2
                    DeviceDiscreMethod = 'Hybrid Trapezoidal';
                case 3
                    DeviceDiscreMethod = 'Virtual Damping';
                otherwise
                    error(['Error: Wrong discretization method.'])
            end
            set_param(gcb,'DiscreMethod',DeviceDiscreMethod);
            set_param(gcb,'LinearizationTimes',num2str(LinearizationTimes));
            if DirectFeedthrough == 1
                set_param(gcb,'DirectFeedthrough','on');
            else
                set_param(gcb,'DirectFeedthrough','off');
            end
            set_param(gcb,'EnableInsideModification','on');
            if DiscreDampingFlag == 1
                set_param(gcb,'DiscreDampingFlag','on');
                set_param(gcb,'DiscreDampingValue',['DeviceDiscreDamping{' num2str(i) '}']);
            else
                set_param(gcb,'DiscreDampingFlag','off');
            end
            set_param(gcb,'EnableInsideModification','off');
            
        end
        
        % If the device is an "infinite bus"
        if floor(DeviceType{i}/10) == 9
            set_param(gcb,'theta0',[num2str(PowerFlow{i}(4))]);
            set_param(gcb,'w','Wbase');
        end
        
    end
end

end