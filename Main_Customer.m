% This function is for customer use

%% Notes

% Please read "README.md" first before using the toolbox.

% For changing custormized data, please use "netlist.xlsx".

% The toolbox defaultly print the results in matlab command window and save
% the results into matlab workspace. Additionally, three figures are
% plotted defaulty: pole map; admittance spectrum measured at each bus;
% frequency related transfer function for each machine.

%% Call "Main"
Main();             % This function runs the toolbox.

%% Custormer available results
GsysDSS;            % Whole-system model (descriptor state space form)
                    % Note: The elements of state, input, and output
                    % vectors are printed in the command window.
                    
GminSS;             % Whole-system model (state space form)
                    % Note: This model is the minimum realization of
                    % GsysDSS, which keeps the same input and output as
                    % GsysDSS, but reduces the order of state.
                    
ListPowerFlow;      % Power flow result
                    % Note: The result is in the form of
                    % | bus | P | Q | V | angle | omega |

ListPowerFlow_;     % Power flow result only for active device by combing 
                    % the load into the nodal admittance matrix.
                    
pole_sys;           % Whole-system poles, or equivalently eigenvalues.

%% Customer function and plot
% Custormer can write their functions here to further deal with the data
% mentioned above.