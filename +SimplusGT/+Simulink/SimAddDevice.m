% This function adds apparatuses into simulink model.

% Author(s): Yitong Li

function [FullName_Device,Name_Device,Pos_Device] = SimAddDevice(Name_Model,Name_LibFile,Size_Device,Shift_Device,Pos_Bus,DeviceBus,DeviceType,ListAdvance)

% Organize data
DiscreMethod = ListAdvance(1);
LinearizationTimes = ListAdvance(2);
DiscreDampingFlag = ListAdvance(3);
DirectFeedthrough = ListAdvance(4);

N_Device = length(DeviceType);

% Add apparatus
for i = 1:N_Device
    if DeviceType{i}~=0100 && DeviceType{i}~=1100
        
        % Get the bus index of apparatus
        Bus = DeviceBus{i};
        
        switch floor(DeviceType{i}/10)
            % ### Ac apparatus
            case 000
                Name_Device{i} = ['SM' num2str(Bus)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Synchronous Machine (dq-Frame System Object)'],FullName_Device{i});
            case 001
                Name_Device{i} = ['VSI-PLL' num2str(Bus)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                if DeviceType{i}~=19
                    add_block([Name_LibFile '/Grid-Following Voltage-Source Inverter (dq-Frame System Object)'],FullName_Device{i});
                else
                    add_block([Name_LibFile '/Grid-Following Inverter (alpha/beta System Object)'],FullName_Device{i});
                end
            case 002
                Name_Device{i} = ['VSI-Droop' num2str(Bus)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Grid-Forming Voltage-Source Inverter (dq-Frame System Object)'],FullName_Device{i});
            case 009
            	Name_Device{i} = ['Inf-Bus' num2str(Bus)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/AC Infinite Bus (Voltage Type)'],FullName_Device{i});
                
            % ### Dc apparatus
            case 101
                Name_Device{i} = ['Buck' num2str(Bus)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Grid-Feeding Buck (System Object)'],FullName_Device{i});
        	case 109
            	Name_Device{i} = ['Inf-Bus' num2str(Bus)];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/DC Infinite Bus (Voltage Type)'],FullName_Device{i});
                
            % ### Interlink
            case 200
              	Name_Device{i} = ['Interlink' num2str(Bus(1)) '-' num2str(Bus(2))];
                FullName_Device{i} = [Name_Model '/' Name_Device{i}];
                add_block([Name_LibFile '/Interlink AC-DC (System Object)'],FullName_Device{i});
                
          	% ### Error check
            otherwise
                error(['Error: DeviceType ' num2str(DeviceType{i}) '.']);
        end
        
        % Set position
       	% The position of apparatus is set by referring to the position of correpsonding bus
        if 0<=DeviceType{i} && DeviceType{i}<=90    
            % For ac apparatuses
            Pos_Device{i} = Pos_Bus{Bus} + Shift_Device;
            set_param(FullName_Device{i},'position',[Pos_Device{i},Pos_Device{i}+Size_Device]);
        elseif 1000<=DeviceType{i} && DeviceType{i}<=1090   
            % For dc apparatuses: smaller
            Pos_Device{i} = Pos_Bus{Bus} + Shift_Device;
            set_param(FullName_Device{i},'position',[Pos_Device{i},Pos_Device{i}+Size_Device-[0,20]]);
        elseif 2000<=DeviceType{i} && DeviceType{i}<=2090
            % For interlink apparatuses: larger
            Pos_Device{i} = Pos_Bus{Bus(1)} + Shift_Device;
            set_param(FullName_Device{i},'position',[Pos_Device{i},Pos_Device{i}+Size_Device+[0,40]]);
        end
        set_param(FullName_Device{i},'Orientation','left');
        
        % Set common variables
      	set_param(gcb,'Sbase','Sbase');
        set_param(gcb,'Vbase','Vbase');
        set_param(gcb,'Ts','Ts');
        
        % For ac apparatus only
        if (DeviceType{i}<1000) || (2000<=DeviceType{i} && DeviceType{i}<=2090)
            set_param(gcb,'Wbase','Wbase');
        end
        
        % For active apparatus only
        if (0<=DeviceType{i} && DeviceType{i}<90) || ...
           (1000<=DeviceType{i} && DeviceType{i}<1090) || ...
           (2000<=DeviceType{i} && DeviceType{i}<2090)
            
            % Set system object parameters
            set_param(gcb,'DeviceType',['DeviceType{' num2str(i) '}']);
            set_param(gcb,'DevicePara',['DevicePara{' num2str(i) '}']);
            set_param(gcb,'PowerFlow',['DevicePowerFlow{' num2str(i) '}']);
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
        
        % If the apparatus is an infinite bus
        if DeviceType{i} == 0090        % Ac infinite bus
            set_param(gcb,'vd',['PowerFlow{' num2str(i) '}(3)']);
            set_param(gcb,'theta0',['PowerFlow{' num2str(i) '}(4)']);
            set_param(gcb,'w','Wbase');
        elseif DeviceType{i} == 1090    % Dc infinite bus
            set_param(gcb,'v',['PowerFlow{' num2str(i) '}(3)']);
        end
        
    end
end

end