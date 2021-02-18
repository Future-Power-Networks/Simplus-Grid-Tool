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
The first type of device is simply a conventional Simulink block. Users just need to copy their custom device into "SimplexPS.xlsx". Then, this device can be found in "Simplex Power System" in the Simulink Library Broswer. This procedure is same to the procedure of [creating custom libraries in Simulink](https://uk.mathworks.com/help/simulink/ug/adding-libraries-to-the-library-browser.html).

### Matlab System Object    
The second type of device is called [Matlab system object](https://uk.mathworks.com/help/simulink/ug/system-design-in-simulink-using-system-objects.html). It defines its input and output relationship in a script, and then links it to a Matlab system block in simulink. This type of device has following features:
 
a) Compared with the Simulink block that can be used in simulations only, Matlab system object can be used by both script (theoratical analysis, e.g., small-signal analysis) and simulink (simulations). This explains why the toolbox can automatically generate the theoratical model and simulation model at the same time.   

b) In simulink, this type of device can run very fast because it supports code generation (such as C).

Adding new device in this form is a bit complicated. The following files need to be changed:

ModelTemplate.m:  
This is the template for adding a new device class, which is the core file. "Inductor.m" and "SynchronousMachine.m" are two specific examples. "Inductor.m" is a very simple example, i.e., a single-phase inductor, which will not be used in the toolbox. "SynchronousMachine.m" is a more complex and practical example, i.e., a three-phase synchronous machine, which is used in toolbox.

SimplexPS.slx:  
This is the simulink library of this toolbox. The simulink block of new device, which is based on the [Matlab System Block](https://uk.mathworks.com/help/simulink/slref/matlabsystem.html), should also be added into the library.

UserData.xlsx:  
The data from users is saved in this file. This file also guides the users about how to use this toolbox. Hence, the introduction of the new devices should be added in this file.

RearrangeDeviceData.m:  
This file gives the default paramters for a new device, and defines how the data in "UserData.xlsx" is transferred to the codes.

DeviceModelCreate.m:  
This file defines the device data transferred from the excel to the class.

SimAddDevice.m:  
This file defines how your device is automatically added when a simulink model is created.