% This function adds transformers into the simulink model.

% Author(s): Yitong Li

function [Name_Trans] = Sim_AddTransformer(Name_Model,Size_Trans,Shift_Trans,Pos_Bus,ListLine,ListSimulation,Shift_ToBus)

% Organize data
fb = ListLine(:,1);     % From bus
tb = ListLine(:,2);     % To bus
Tbr  = ListLine(:,7);   % Transformer tap ratio
N_Branch = length(fb);

F0 = ListSimulation(length(ListSimulation));

Name_Trans{1} = [];

% Add transformer
for i = 1:N_Branch
    if fb(i)~= tb(i)
        if Tbr(i) ~= 1
            Pos_Trans{i} = [Pos_Bus{tb(i)}(1),Pos_Bus{fb(i)}(2)];

            Name_Trans{i} = ['Transformer' num2str(fb(i)) num2str(tb(i))];
            FullName_Trans{i} = [Name_Model '/' Name_Trans{i}];
            add_block(['powerlib/Elements/Three-Phase Transformer (Two Windings)'],FullName_Trans{i});

            % Set position
            Pos_Trans{i} = Pos_Trans{i} + [Shift_Trans(1)*Shift_ToBus{i},Shift_Trans(2)];
            set_param(gcb,'position',[Pos_Trans{i},Pos_Trans{i}+Size_Trans]);
            set_param(gcb,'Orientation','down');
            set_param(gcb,'Measurements','None');

            % Set parameters
            set_param(gcb,'NominalPower',['[1, ' num2str(F0) ']']);
            set_param(gcb,'Winding1',['[' num2str(Tbr(i)) ', 1e-9, 1e-9]']);
            set_param(gcb,'Winding2','[1, 1e-9, 1e-9]');
            set_param(gcb,'Rm','1e4');
            set_param(gcb,'Lm','1e4');
        end
    end
end

end