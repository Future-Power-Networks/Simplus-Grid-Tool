Here are the two models in our paper. I am using matlab 2018b. But other versions should also be fine.


1) Three-node circuit:
 add the foler '...\three_node_circuit' into the path of matlab, then run 'Case_SimpleCircuit_3nodes.m', results are saved in the variables in the workspace.

2) IEEE 14-buS system:
The original data comes from the dynamic model of University of Cyprus (KIOS research team), and the original model is built in DigSLIENT PowerFactory. Here I transplanted the model in Simplus Grid-tool, and the two are compared to be nearly the same on small-signal dynamics (the pole-map are compared to be the same). Then three additional IBRs are added in our model to verify sensitivity method in our paper.
Instructions to run the model:
Add '...\modified_ieee_14\' folder into the path, including all subfolders.
Then run 'UserMain.m'. You can read those comments in the code if you wish to acuqire more information.
If you are going to build your own model using Simplus Grid-tool, it is recommonded that you fork from 'master' branch of Simplus Grid-tool instead of this branch.

Contact: yue.zhu18@imperial.ac.uk