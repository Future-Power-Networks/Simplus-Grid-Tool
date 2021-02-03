# Advanced Manual for Developers

============  
Introduction  
============  

This is the manual will introduce more details of developing this toolbox. Please read "DeveloperManual.md" first before reading this one.

============  
Adding New Device  
============

Files need to be changed for adding a new device:

UserMain.xlsx:  
The data from users is saved in this file. This file also guides the users about how to use this toolbox. The introduction about the new devices should be added into this file.

ModelTemplate.m:  
This is the template for adding a new device class.

RearrangeDeviceData.m:  
This file defines the default paramters for a new device.

DeviceModelCreate.m:  
This file defines the device data transferred from the excel to the class.

SimplexPS.slx:  
This is the library of this toolbox. The simulink block of new device should also be added into the library.

SimAddDevice.m:  
This file defines how your device is added automatically when a simulink model is created.