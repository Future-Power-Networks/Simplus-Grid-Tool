% This function adds bus into simulink model

% Author(s): Yitong Li

function [Name_Bus,Pos_Bus]= SimAddBus(Name_Model,Name_LibFile,Size_Bus,Pos_Bus,ListBus,Dist_Bus)

% Organize data
N_Bus = length(ListBus(:,1));

% Add blocks
for i = 1:N_Bus
    
    [~,~,AreaType{i}] = SimplusGT.Toolbox.CheckBus(i,ListBus);
    
    Name_Bus{i} = ['Bus' num2str(i)];
    FullName_Bus{i} = [Name_Model '/' Name_Bus{i}];
    if AreaType{i} == 1
        add_block([Name_LibFile '/Three-Phase Bus'],FullName_Bus{i});   % Add ac bus
    elseif AreaType{i} == 2
        add_block([Name_LibFile '/Two-Phase Bus'],FullName_Bus{i});     % Add dc bus
    end
    set_param(gcb,'position',[Pos_Bus{i},Pos_Bus{i}+Size_Bus]);
    Pos_Bus{i+1} = Pos_Bus{i} + Dist_Bus;

end

end