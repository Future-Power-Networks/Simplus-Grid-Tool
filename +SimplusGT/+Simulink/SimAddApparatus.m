% This function adds apparatuses into simulink model.

% Author(s): Yitong Li

function [FullName_Apparatus,Name_Apparatus,Pos_Apparatus] = SimAddApparatus(Name_Model,Name_LibFile,Size_Apparatus,Shift_Apparatus,Pos_Bus,ApparatusBus,ApparatusType,Advance)

% Organize data
DiscreMethod = Advance.DiscretizationMethod;
LinearizationTimes = Advance.LinearizationTimes;
DiscreDampingFlag = Advance.DiscretizationDampingFlag;
DirectFeedthrough = Advance.DirectFeedthrough;

N_Apparatus = length(ApparatusType);

% Add apparatus
for i = 1:N_Apparatus
    if ApparatusType{i}~=0100 && ApparatusType{i}~=1100
        
        % Get the bus index of apparatus
        Bus = ApparatusBus{i};
        
        switch floor(ApparatusType{i}/10)
            % ### Ac apparatus
            case 000
                Name_Apparatus{i} = ['SM' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Synchronous Machine (dq-Frame System Object)'],FullName_Apparatus{i});
            case 001
                Name_Apparatus{i} = ['VSI-PLL' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                if ApparatusType{i}~=19
                    add_block([Name_LibFile '/Grid-Following Voltage-Source Inverter (dq-Frame System Object)'],FullName_Apparatus{i});
                else
                    add_block([Name_LibFile '/Grid-Following Inverter (alpha/beta System Object)'],FullName_Apparatus{i});
                end
            case 002
                Name_Apparatus{i} = ['VSI-Droop' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Grid-Forming Voltage-Source Inverter (dq-Frame System Object)'],FullName_Apparatus{i});
            case 003
                Name_Apparatus{i} = ['BESS' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Battery (dq-Frame System Object)'],FullName_Apparatus{i});
            case 004
                if ApparatusType{i}==40
                    Name_Apparatus{i} = ['PV_GFM' num2str(Bus)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/Photovoltaic-GFM (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i}==41
                    Name_Apparatus{i} = ['PV_GFL' num2str(Bus)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/Photovoltaic-GFL (dq-Frame System Object)'],FullName_Apparatus{i});
                end
            case 005
                if ApparatusType{i} == 50
                    Name_Apparatus{i} = ['Wind-GFM' num2str(Bus)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/Grid-Forming Wind Turbine Generator System (dq-Frame System Object)'],FullName_Apparatus{i});
                elseif ApparatusType{i} == 51
                    Name_Apparatus{i} = ['Wind-GFL' num2str(Bus)];
                    FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                    add_block([Name_LibFile '/Grid-Following Wind Turbine Generator System (dq-Frame System Object)'],FullName_Apparatus{i});
                end
            case 009
            	Name_Apparatus{i} = ['Inf-Bus' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/AC Infinite Bus (Voltage Type)'],FullName_Apparatus{i});
                
            % ### Dc apparatus
            case 101
                Name_Apparatus{i} = ['Buck' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Grid-Feeding Buck (System Object)'],FullName_Apparatus{i});
        	case 109
            	Name_Apparatus{i} = ['Inf-Bus' num2str(Bus)];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/DC Infinite Bus (Voltage Type)'],FullName_Apparatus{i});
                
            % ### Interlink
            case 200
              	Name_Apparatus{i} = ['Interlink' num2str(Bus(1)) '-' num2str(Bus(2))];
                FullName_Apparatus{i} = [Name_Model '/' Name_Apparatus{i}];
                add_block([Name_LibFile '/Interlink AC-DC (System Object)'],FullName_Apparatus{i});
                
          	% ### Error check
            otherwise
                error(['Error: ApparatusType ' num2str(ApparatusType{i}) '.']);
        end
        
        % Set position
       	% The position of apparatus is set by referring to the position of correpsonding bus
        if 0<=ApparatusType{i} && ApparatusType{i}<=90    
            % For ac apparatuses
            Pos_Apparatus{i} = Pos_Bus{Bus} + Shift_Apparatus;
            set_param(FullName_Apparatus{i},'position',[Pos_Apparatus{i},Pos_Apparatus{i}+Size_Apparatus]);
        elseif 1000<=ApparatusType{i} && ApparatusType{i}<=1090   
            % For dc apparatuses: smaller
            Pos_Apparatus{i} = Pos_Bus{Bus} + Shift_Apparatus;
            set_param(FullName_Apparatus{i},'position',[Pos_Apparatus{i},Pos_Apparatus{i}+Size_Apparatus-[0,20]]);
        elseif 2000<=ApparatusType{i} && ApparatusType{i}<=2090
            % For interlink apparatuses: larger
            Pos_Apparatus{i} = Pos_Bus{Bus(1)} + Shift_Apparatus;
            set_param(FullName_Apparatus{i},'position',[Pos_Apparatus{i},Pos_Apparatus{i}+Size_Apparatus+[0,40]]);
        end
        set_param(FullName_Apparatus{i},'Orientation','left');
        
        % Set common variables
      	set_param(gcb,'Sbase','Sbase');
        set_param(gcb,'Vbase','Vbase');
        set_param(gcb,'Ts','Ts');
        
        % For ac apparatus only
        if (ApparatusType{i}<1000) || (2000<=ApparatusType{i} && ApparatusType{i}<=2090)
            set_param(gcb,'Wbase','Wbase');
        end
        
        % For active apparatus only
        if (0<=ApparatusType{i} && ApparatusType{i}<90) || ...
           (1000<=ApparatusType{i} && ApparatusType{i}<1090) || ...
           (2000<=ApparatusType{i} && ApparatusType{i}<2090)
            
            % Set system object parameters
            set_param(gcb,'ApparatusType',['ApparatusType{' num2str(i) '}']);
            set_param(gcb,'ApparatusPara',['ApparatusPara{' num2str(i) '}']);
            set_param(gcb,'PowerFlow',['ApparatusPowerFlow{' num2str(i) '}']);
            set_param(gcb,'x0',['x_e{' num2str(i) '}']);
            set_param(gcb,'OtherInputs',['OtherInputs{' num2str(i) '}']);

            % Set discretization method
            switch DiscreMethod
                case 1
                    ApparatusDiscreMethod = 'Forward Euler';
                case 2
                    ApparatusDiscreMethod = 'Hybrid Trapezoidal';
                case 3
                    ApparatusDiscreMethod = 'Virtual Damping';
                otherwise
                    error(['Error: Wrong discretization method.'])
            end
            set_param(gcb,'DiscreMethod',ApparatusDiscreMethod);
            set_param(gcb,'LinearizationTimes',num2str(LinearizationTimes));
            if DirectFeedthrough == 1
                set_param(gcb,'DirectFeedthrough','on');
            else
                set_param(gcb,'DirectFeedthrough','off');
            end
            set_param(gcb,'EnableInsideModification','on');
            if DiscreDampingFlag == 1
                set_param(gcb,'DiscreDampingFlag','on');
                set_param(gcb,'DiscreDampingValue',['ApparatusDiscreDamping{' num2str(i) '}']);
            else
                set_param(gcb,'DiscreDampingFlag','off');
            end
            set_param(gcb,'EnableInsideModification','off');
            
        end
        
        % If the apparatus is an infinite bus
        if ApparatusType{i} == 0090        % Ac infinite bus
            set_param(gcb,'vd',['PowerFlow{' num2str(i) '}(3)']);
            set_param(gcb,'theta0',['PowerFlow{' num2str(i) '}(4)']);
            set_param(gcb,'w','Wbase');
        elseif ApparatusType{i} == 1090    % Dc infinite bus
            set_param(gcb,'v',['PowerFlow{' num2str(i) '}(3)']);
        end
        
    end
end

end