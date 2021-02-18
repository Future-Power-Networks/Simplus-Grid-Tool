# Manual for Developers

============  
Introduction  
============

This is the quick manual for developers who want to contribute to the toolbox developement. If you have any inquires, please feel free to contact leading developers: Yitong Li (yitong.li15@imperial.ac.uk), Yunjie Gu (yg934@bath.ac.uk), and Yue Zhu (yue.zhu18@imperial.ac.uk).

============  
Tips for Developers  
============  

1- The toolbox development is managed by GitHub. For getting involved, the developers should register a GitHub account first and [fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) and star this toolbox respository. [GitHub Desktop](https://desktop.github.com/) is also recommended, which helps to manage/update the respository. Only the leading developers have the access to this toolbox master respository. But everyone can change their forked respository (which is essentially a cloned version under your own account), and can raise a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests). Leading developers will then receive a notification for reviewing the modified codes in the pull request, and decide to accept it (i.e., merge your modified codes to the master respository) or reject it. This is the formal way of collaboration in GitHub.

2- The toolbox structure will be introduced briefly below. After that, you can also view the advanced manual in "Files for Developers" subfolder, and/or view the codes and comments in each ".m" file for more details.

3- There are three sets of naming conventions:  
    1) microsoft, like PlotFigure, EpwmConfig  
    2) gnu-linux, like plot_figure, pwm_config  
    3) free style, like ss2tf  
Matlab uses free style for its common functions, and uses microsoft mainly in others. We also use microsoft mainly in this project, and use free style when needed.

4- The project is going to be renamed SimplexPowerSystem (SimplexPS, in short), while the project on HIL is called SimplexRealTime (SimplexRT).

============  
File structure  
============  

The file structure (documentation tree) of the toolbox is shown below. The detailed introduction of each folder will be given next.

<pre>
Root path  
	-> +SimplexPS (root name space)  
		-> +Class (sub name sapce, class of matlab system object)  
		-> +PowerFlow (sub name space, power flow analysis)  
		-> +Toolbox (sub name space, toolbox basic functions)  
		-> +Simulink (sub name space, simulink functions)  
		-> Generic functions (such as "dss2ss".)  
	-> Debug (functions used by developers for debug.)  
	-> Documentations (manuals for users and developers)  
	-> Examples  
	-> Library (library for Simulink blocks)  
		-> slblocks.m (load SimplexPS lib to Simulink broswer)  
		-> SimplexPS.slx (library file)  
		-> SimplexPS_2015a.xls (library file for matlab 2015a)  
	-> InstallSimplexPS.m  
	-> UninstallSimplexPS.m  
	-> UserMain.m  
	-> UserData.xlsx  
	-> README.md
</pre>

============  
Root Path  
============  

Root path contains the main folders for toolbox and the files for users.

### Files for both users and developers:

"README.m" is a quick introduction of the toolbox. It should be changed by leading developers only.

"InstallSimplexPS.m" and "UninstallSimplexPS.m" are files for the installation and uninstallation of the toolbox.

"UserMain.m" is the main file for users, which runs the toolbox. "UserDate.xlsx" is file containing all data of the users' power systems. Default data of a 4-bus power system is saved in this excel file, so that users can run "UserMain.m" directly the first time without the requirement of any further actions or inputting their own data.

"Examples" folder contains the examples which can be used by users. For example, "IEEE_14Bus.xlsx" contains the data for a standard IEEE 14 bus power system.

### Folders for developers only:

"+SimplexPS" is the root name space folder, which contains all the functions (except for functions available for users such as "UserMain.m") used in the toolbox. For Matlab, "+" implies [a name space (or equivalently a package)](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html). The use of name space is mainly for avoiding the unexpected duplication of function names and errors. Functions are further classfied into different groups and saved in different sub name spaces, as introduced in detail later.

"Debug" folder contains the backup files used by different developers for debugging the toolbox. These files are for debug only and should NOT be used/called when users run the toolbox.

"Documentations" folder contains documentations which introduce the toolbox, such as "DeveloperManual.md" and "UserManual.md".

"Library" folder contains the library files for Simulink blocks used by the toolbox. "SimplexPS_2015a.xls" is the library file in Matlab version 2015a. Runing "InstallSimplexPS.m" will automatically convert this lib file to the version of users' Maltab.

============  
Root Name Space: "+SimplexPS" Folder  
============  

This is the root name space folder of the toolbox, and contains all functions used by the toolbox (except for functions available for users directly such as "UserMain.m"). Generic functions (used for generic purpose such as mathmatical calculations, bode plot, etc) are saved in this root name space folder directly. Other advanced functions are saved in corresponding sub name space folders ("+Toolbox", "+Simulink", etc), as introduced next.

### "+Toolbox" Sub-folder

This folder contains the basic functions for toolbox. For example, "Main.m" is the main function for toolbox; "DeviceModelCreate.m" creates the small-signal state space models of devices.

### “+Simulink” Sub-folder

This folder contains the functions generating the Simulink model from user data programmatically and automatically. "MainSimulinkModel.m" is the main function for this folder.

### "+Class" Sub-folder

This folder contains the class functions of Matlab system object, which defines the models of devices used for both small-signal modeling in Matlab script and time-domain simulation in Simulink.

The class tree is shown here:

<pre>
matlab.System  
	-> ModelBase  
		& matlab.system.mixin.Nondirect  
		& matlab.system.mixin.Propagates  
		-> ModelAdvance  
			-> ModelTemplate  
			-> SynchronousMachine  
			-> GridFollowingVSI  
			-> GridFormingVSI  
			-> InfiniteBus  
			-> FloatingBus  
			-> Inductor  
</pre>

"matlab.System" is the base class of Matlab itself. "matlab.system.mixin.Nondirect" and "matlab.system.mixin.Propagates" are also classes of Matlab, which are included for using system objects in Simulink.

"ModelBase.m" is the base class of the toolbox, which defines the base properties (such as state space matrics, etc) and base methods (such as the linearization algorithms, the functions for reading/writing properties, etc).

"ModelAdvance.m" is the sub-class of "ModelBase.m". It is the advanced class used for Matlab system block in Simulink. It defines advanced properties and methods (such as discretization methods).

"SynchronousMachine.m" is the specific device class for synchronous machines, which contains the detailed physical models of synchronous machines (such as signal lists, equilibrium calculations, state equations, etc). Similarly, "GridFollowingVSI.m" is the specific device class for grid following voltage source inverters ...

"ModelTemplate.m" is the template class for a specific device. Developers can use this file as a starting point to build/customize a new device.

### "+PowerFlow" Sub-folder

This folder contains functions for power flow analysis. "PowerFlowGS.m" and "PowerFlowNR" correspond to Gauss-Seidel and Newton–Raphson methods respectively.

### "+Modal" Sub-folder

This folder contains functions for modal analysis (such as participation analysis, sensitivity analysis).