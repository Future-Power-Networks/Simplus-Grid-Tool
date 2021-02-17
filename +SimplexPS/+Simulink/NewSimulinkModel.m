% This functon creates a new simulink model

% Author(s): Yitong Li

%% Notes
%
% Can use "get_param" to check the properties.

%%
function NewSimulinkModel(varargin)

%% Get the parameters
for n = 1:length(varargin)
    if(strcmpi(varargin{n},'ModelName'))
        ModelName = varargin{n+1};
     elseif(strcmpi(varargin{n},'Solver'))
         Solver = varargin{n+1};
    elseif(strcmpi(varargin{n},'ScreenColor'))
        ScreenColor = varargin{n+1};
    end
end

%% Set the default values
try   
    ModelName;
catch  
    ModelName = 'my_untitled';
end

try
    Solver;
catch
    Solver = 'FixedStepAuto';
end

try
    FixedStep;
catch
    FixedStep = 'Ts';
end

try
    StartTime;
catch
    StartTime = '0';
end

try
    StopTime;
catch
    StopTime = '5';
end

try
    ScreenColor;
catch
    ScreenColor = 'white';
end

%% Create a blank model
% Create a model
t2 = new_system(ModelName);

% Open the model
open_system(ModelName);

% Set screen color
set_param(ModelName,'ScreenColor',ScreenColor);

% Set solver
set_param(ModelName,'Solver',Solver);
set_param(ModelName,'FixedStep',FixedStep);

% Set start time and stop time
set_param(ModelName,'StartTime',StartTime);
set_param(ModelName,'StopTime',StopTime);

%% Save the model
% save_system(ModelName);

end