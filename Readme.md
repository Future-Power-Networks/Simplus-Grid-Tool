Here are the two case studies in our paper. I am using Matlab 2018b. But other versions should also be fine.

* Three-node circuit:
 add the foler 'three_node_circuit' into the path of Matlab, then run 'Case_SimpleCircuit_3nodes.m', results are saved in the variables in the workspace.

* IEEE 14-buS system:
The original dynamic model comes from [here](https://www2.kios.ucy.ac.cy/testsystems/index.php/ieee-14-bus-modified-test-system/), which is built in DigSLIENT PowerFactory. Here I transplanted the model into Simplus Grid-tool, and the two are compared to be nearly the same for small-signal dynamics (the pole-maps are the same). Three additional IBRs are further added in my model to verify the sensitivity method in our paper.
Instructions:
add the folder 'modified_ieee_14' folder into the path, including all subfolders, then run 'UserMain.m'. You can read those comments in the code if you wish to acuqire more information. If you are going to build your own model using Simplus Grid-tool, it is recommonded that you fork from 'master' branch of this repo.

Contact: yue.zhu18@imperial.ac.uk
