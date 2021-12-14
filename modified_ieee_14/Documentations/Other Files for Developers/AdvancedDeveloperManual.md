# Advanced Manual for Developers

## Introduction

This manual will introduce more details of developing this toolbox. Please read "DeveloperManual.md" first before reading this one.

## Adding New Apparatus

There are two different types of apparatuses in the toolbox:    
1- Simulink block;    
2- Matlab system object.

### Simulink Block   
The first type of apparatus is simply a conventional Simulink block. Users just need to copy their custom apparatus into "SimplusGT.slx". Then, this apparatus can be found in "Simplus Grid Tool" in the Simulink Library Broswer. This procedure is same to the procedure of [creating custom libraries in Simulink](https://uk.mathworks.com/help/simulink/ug/adding-libraries-to-the-library-browser.html).

### Matlab System Object    
The second type of apparatus is called [Matlab system object](https://uk.mathworks.com/help/simulink/ug/system-design-in-simulink-using-system-objects.html). It defines its input and output relationship in a script, and then links it to a Matlab system block in simulink. This type of apparatus has following features:
 
a) Compared with the Simulink block that can be used in simulations only, Matlab system object can be used by both script (theoratical analysis, e.g., small-signal analysis) and simulink (simulations). This explains why the toolbox can automatically generate the theoratical model and simulation model at the same time.   

b) In simulink, this type of apparatus can run very fast because it supports code generation (such as C).

Adding new apparatus in this form is a bit complicated. The following files need to be changed:

ModelTemplate.m:  
This is the template for adding a new apparatus class, which is the core file. "Inductor.m" and "SynchronousMachine.m" are two specific examples. "Inductor.m" is a very simple example, i.e., a single-phase inductor, which will not be used in the toolbox. "SynchronousMachine.m" is a more complex and practical example, i.e., a three-phase synchronous machine, which is used in toolbox.

SimplusGT.slx:  
This is the simulink library of this toolbox. The simulink block of new apparatus, which is based on the [Matlab System Block](https://uk.mathworks.com/help/simulink/slref/matlabsystem.html), should also be added into the library.

UserData.xlsx:  
The data from users is saved in this file. This file also guides the users about how to use this toolbox. Hence, the introduction of the new apparatuses should be added in this file.

RearrangeApparatusData.m:  
This file gives the default paramters for a new apparatus, and defines how the data in "UserData.xlsx" is transferred to the codes.

ApparatusModelCreate.m:  
This file defines the apparatus data transferred from the excel to the class.

SimAddApparatus.m:  
This file defines how your apparatus is automatically added when a simulink model is created.