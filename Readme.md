Here are the three case studies in our paper. I am using Matlab 2018b. But other versions should also be fine.

* Three-node circuit:
 add the foler 'three_node_circuit' into the path of Matlab, then run 'Case_SimpleCircuit_3nodes.m', results are saved in the variables in the workspace.

* IEEE 14-bus system:
The original dynamic model comes from [here](https://www2.kios.ucy.ac.cy/testsystems/index.php/ieee-14-bus-modified-test-system/), which is built in DigSLIENT PowerFactory. Here I transplanted the model into Simplus Grid-tool, and the two are compared to be nearly the same for small-signal dynamics (the pole-maps are the same). Three additional IBRs are further added in my model to verify the sensitivity method in our paper.
Instructions:
add the folder 'system_cases' folder into the path, including all subfolders, then in 'UserMain.m' select the case index and run the script. You can read those comments in the scripts for more information. If you are going to build your own model using Simplus Grid-tool, it is recommonded that you fork from 'master' branch of this repo.

* NETS-NYPS 68-bus system:
add the folder 'system_cases' folder into the path, including all subfolders, then in 'UserMain.m' select the case index and run the script.


Contact: yue.zhu18@imperial.ac.uk
