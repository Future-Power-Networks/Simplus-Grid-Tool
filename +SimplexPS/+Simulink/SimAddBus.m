% This function adds bus into simulink model

% Author(s): Yitong Li

function [Name_Bus,Pos_Bus]= SimAddBus(Name_Model,Name_LibFile,Size_Bus,Pos_Bus,ListLine,Dist_Bus)

% Organize data
fb = ListLine(:,1); % From bus
tb = ListLine(:,2); % To bus
N_Bus = max(max(fb),max(tb));

% Add blocks
for i = 1:N_Bus
    
    Name_Bus{i} = ['Bus' num2str(i)];
    FullName_Bus{i} = [Name_Model '/' Name_Bus{i}];
    add_block([Name_LibFile '/Three-Phase Bus'],FullName_Bus{i});
    set_param(gcb,'position',[Pos_Bus{i},Pos_Bus{i}+Size_Bus]);
    Pos_Bus{i+1} = Pos_Bus{i} + Dist_Bus;

end

end