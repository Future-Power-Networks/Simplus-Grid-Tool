function GreyboxExcelGenerate(N_Device,DeviceType,DeviceStateStr,DeviceInputStr,...
DeviceOutputStr,ZbusStateStr,GsysDSS)

xlswrite('GreyBoxConfig.xlsx',DeviceStateStr{1},'state','B2');

end