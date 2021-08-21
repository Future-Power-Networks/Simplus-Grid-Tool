% find variables from varargin

function [NewVar,flag] = LoadVar(Var,Name,varargin_)

flag = 0;
for n = 1:length(varargin_)
    if iscell(varargin_{n})
        [NewVar,flag] = LoadVar(Var,Name,varargin_{n});
        if flag == 1
            break;
        end
    else
        %if(strcmpi(varargin_{n},Name))
        if(strcmp(varargin_{n},Name))
            NewVar = varargin_{n+1};
            flag = 1;
            break;
        end
    end
end

if flag == 0
    NewVar = Var;
end

end