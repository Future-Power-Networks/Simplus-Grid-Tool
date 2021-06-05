% This function adds transformers into the simulink model.

% Author(s): Yitong Li

function [Name_Trans] = SimAddTransformer(Name_Model,Size_Trans,Shift_Trans,Pos_Bus,ListLine,Shift_ToBus)

% Organize data
FB = ListLine(:,1);         % From bus
TB = ListLine(:,2);         % To bus
Tbr  = ListLine(:,7);       % Transformer tap ratio
AreaType = ListLine(:,8);
N_Branch = length(FB);

Name_Trans{1} = [];

% Add transformer
for i = 1:N_Branch
    if AreaType == 1        % Add transformer for ac branch only
        if FB(i)~= TB(i)
            if Tbr(i) ~= 1
                Pos_Trans{i} = [Pos_Bus{TB(i)}(1),Pos_Bus{FB(i)}(2)];

                % Get the block name
                Name_Trans{i} = ['Transformer' num2str(FB(i)) '-' num2str(TB(i))];
                for j = 1:(length(Name_Trans)-1)
                    if FB(i)==FB(j) && TB(i)==TB(j)
                        Name_Trans{i} = [Name_Trans{i},'_'];    % Avoid the name duplication
                    end
                end
                FullName_Trans{i} = [Name_Model '/' Name_Trans{i}];

                % Add block
                add_block(['powerlib/Elements/Three-Phase Transformer (Two Windings)'],FullName_Trans{i});

                % Set position
                Pos_Trans{i} = Pos_Trans{i} + [Shift_Trans(1)*Shift_ToBus{i},Shift_Trans(2)];
                set_param(gcb,'position',[Pos_Trans{i},Pos_Trans{i}+Size_Trans]);
                set_param(gcb,'Orientation','down');
                set_param(gcb,'Measurements','None');

                % Set parameters
                set_param(gcb,'NominalPower','[Sbase, Fbase]');
                set_param(gcb,'Winding1',['[' num2str(Tbr(i)) '*Vbase, 0, 0]']);
                set_param(gcb,'Winding2','[Vbase, 0, 0]');
                set_param(gcb,'Rm','1e4');
                set_param(gcb,'Lm','1e4');
            end
        end
    end
end

end