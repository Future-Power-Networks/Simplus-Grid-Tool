====================================
Document tree:
====================================

Rootpath
	-> +SimplexPS (root name space)
		-> +Class
		-> +PowerFlow
		-> +Toolbox
		-> +Simulink
		-> Generic functions such as dss2ss ...
	-> Documentations
	-> Examples
	-> Library (Simulink library for SimplexPS)
		-> slblocks.m (Load SimplexPS lib to Simulink broswer)
		-> SimplexPS.slx (SimplexPS libaray file)
	-> InstallSimplexPS.m
	-> UnstallSimplexPS.m
	-> CustomerMain.m
	-> CustomerData.xlsx
	-> README.md

====================================
Notes:
====================================

1. The Project is going to be renamed SimplexPowerSystem (SimplexPS, in short), while the project on HIL is called SimplexRealTime (SimplexRT)

2. under the root folder there will be three main folders: +SimplexPS (root name space), Examples, and Documentations

3. +SimplexPS will contain all source codes and sub-namespaces.

4. there are three set of naming conventions: 
    1) microsoft, like PlotFigure EpwmConfig
    2) gnu-linux, like plot_figure pwm_config
    3) free style, like ss2tf
    Matlab use free style for its common functions, and use microsoft in others. Let's use microsoft mainly, and use free style when needed
    e.g. when we override a matlab function, or defining similar general functions, we can use free style. For all others, we use microsoft

5. under ther project root folder, there should be an install.m and uninstall.m, with set the root folder in to matlab path, and unset it, respectively, beside other necessary routines needed to get the enviroment ready for users: say resaving the SimplexPS.slx lib file in the version of the running PC.
