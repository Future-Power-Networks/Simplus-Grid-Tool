# Advanced Manual for Developers

============  
Introduction  
============  

This is the manual will introduce more details of developing this toolbox. Please read "DeveloperManual.md" first before reading this one.

============  
Adding New Device  
============

There are two different types of devices in the toolbox:    
1- Simulink block;    
2- Matlab system object.

### Simulink Block   
The first type of device is simply a conventional Simulink block. Users just need to copy their custom device into "SimplexPS.xlsx". Then, this device can be found in "Simplex Power System" in the Simulink Library Broswer.This procedure is same to the procedure of [creating custom libraries in Simulink](https://uk.mathworks.com/help/simulink/ug/adding-libraries-to-the-library-browser.html).

### Matlab System Object    
The second type of device is called [Matlab system object](https://uk.mathworks.com/help/simulink/ug/system-design-in-simulink-using-system-objects.html).It defines its input and output relationship in a class script, and then links it to a Matlab system block in simulink. This type of device has following features:   
 
a) Compared with the Simulink block that can be used in simulations only, Matlab system object can be used by both script (theoratical analysis, e.g., small-signal analysis) and simulink (simulations). This is the main reason why we use this type of device in the toolbox, so that the toolbox can automatically generate the theoratical model and simulation model at the same time.   

b) In simulink, this type of device can run very fast because it supports code generation (such as C).

Files need to be changed for adding a new device:

UserData.xlsx:  
The data from users is saved in this file. This file also guides the users about how to use this toolbox, for example, the introduction of the new devices should be added into this file.

ModelTemplate.m:  
This is the template for adding a new device class, which is the core file.

RearrangeDeviceData.m:  
This file defines the default paramters for a new device.

DeviceModelCreate.m:  
This file defines the device data transferred from the excel to the class.

SimplexPS.slx:  
This is the simulink library of this toolbox. The simulink block of new device should also be added into the library.

SimAddDevice.m:  
This file defines how your device is added automatically when a simulink model is created.