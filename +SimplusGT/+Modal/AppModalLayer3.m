% This function is for the API of layer 3.
% 
% Author(s): Yue Zhu
% Modified by: Yitong Li

function MdLayer3=AppModalLayer3(AppSelect, MdDataSave)

ApparatusType = evalin('base', 'ApparatusType');
Para = evalin('base', 'Para');
ApparatusPowerFlow = evalin('base', 'ApparatusPowerFlow');
Ts = evalin('base', 'Ts');
ApparatusBus = evalin('base', 'ApparatusBus');
ListBus = evalin('base', 'ListBus');


ResidueAll=MdDataSave.ResidueAll;
ZmValAll = MdDataSave.ZmValAll;
MdMode = MdDataSave.MdMode;
ModeSelAll = MdDataSave.ModeSelAll;

Mode_Hz = MdMode(ModeSelAll);
MdLayer3 = SimplusGT.Modal.MdLayer3(ResidueAll{1},ZmValAll{1},Mode_Hz,ApparatusType,...
        AppSelect,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus);

end
