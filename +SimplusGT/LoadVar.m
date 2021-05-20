% This function find variables from varargin

% Author(s): Yunjie Gu, Yitong Li

%%
% Notes:
%
% Function input:
% Var       - default value
% Name      - string name of the target variable
% varargin_ - varargin input

%%
function [NewVar,Flag] = LoadVar(Var,Name,varargin_)

% Flag is used to check if NewVar is set by varargin or by default Var
Flag = 0;

for n = 1:length(varargin_)
    if iscell(varargin_{n})
        %error(['Error: varargin{n} is cell-type data']);
        [NewVar,flag] = SimplusGT.LoadVar(Var,Name,varargin_{n});
        if flag == 1
            break;
        end
    else
        % Find the target variable
        if(strcmp(varargin_{n},Name))
            NewVar = varargin_{n+1};
            Flag = 1;   % Target variable was found sucessfully
            break;
        end
    end
end

% Set NewVar by a default value
if Flag == 0
    NewVar = Var;
end

end