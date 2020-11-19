% This function adds powergui into the simulink model

% Author(s): Yitong Li

function SimAddPowerGUI(Name_Model,Size_powergui,Position_powergui)

FullName_powergui = [Name_Model '/powergui'];
add_block(['powerlib/powergui'],FullName_powergui);
set_param(gcb,'position',[Position_powergui,Position_powergui+Size_powergui]);
set_param(gcb,'SimulationMode','Discrete');
set_param(gcb,'SampleTime','Ts');

end