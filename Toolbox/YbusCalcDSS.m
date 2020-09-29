% This function calculate the descriptor (implicit) state space model of
% nodal admittance matrix.

% Author(s): Yitong Li, Yunjie Gu

%% Notes

% Self-branch: branches with from = to, connecting the bus and ground
% Mutual-branch: branches with from ~= to, connecting two difference buses

% FrameFlag: 1-abc, 2-alpha/beta, 3-dq

%% Function

function [YbusObj,YbusDSS,YbusCell] = YbusCalcDSS(linedata,w) 

% Load the data
fb = linedata(:,1);             % From bus number...
tb = linedata(:,2);             % To bus number...
Rlist = linedata(:,3);              % Resistance,  R...
Xlist = linedata(:,4);              % Inductance,  wL...
Blist = linedata(:,5);              % Capacitance, wC...
Glist = linedata(:,6);              % Conductance, G...

n_bus = max(max(fb),max(tb));           % Number of buses
n_branch = length(fb);                  % Number of branches, including self branches
      
%%
% Calculate the state space model of each branch
for n = 1:n_branch              % Calculate the branch paramter one by one

    R = Rlist(n);
    X = Xlist(n);
    B = Blist(n);
    G = Glist(n);

    if ( isinf(R) || isinf(X) || ( (G==0) && (B==0) ) )     % open circuit
        A_op = []; B_op = []; E_op = [];
        C_op = []; D_op = [0,0;0,0];
        Ybranch{n} = dss(A_op,B_op,C_op,D_op,E_op); 
    elseif ( (R==0) && (X==0) && ( isinf(G) || isinf(B) ) )
        error(['Error: short circuit, branch from ' num2str(fb(n)) ' to ' num2str(tb(n))]);
    elseif ( isinf(G) || isinf(B) )  	% RL branch, normally for mutual branch
        % KVL equations
        % [vd] = {[R  0] + [sL -wL]}*[id]
        % [vq]    [0  R]   [wL  sL]  [iq]
        if X == 0                       % R branch
            A_RL = []; B_RL = []; E_RL = [];
            C_RL = []; D_RL = inv([R,0;0,R]);
        else                            % RL or L branch
            % State equations
            % [d_id]/dt = 1/L*[-R  wL]*[id] + 1/L*[1 0]*[vd]
            % [d_iq]          [-wL -R] [iq]       [0 1] [vq]
            A_RL = 1/(X/w) * [-R,X;-X,-R];
            B_RL = 1/(X/w) * [1,0;0,1];
            E_RL = [1,0;0,1];
            % Output equations
            % [id] = [1 0]*[id] + [0 0]*[vd]
            % [iq] = [0 1] [iq]   [0 0] [vq]
            C_RL = [1,0;0,1];
            D_RL = [0,0;0,0];
        end
        % Get the branch model
        Y_RL = dss(A_RL,B_RL,C_RL,D_RL,E_RL); 
        Ybranch{n} = Y_RL;
    elseif ( (R==0) && (X==0) )         % GC branch, normally for self branch
        % KCL equations
        % [id] = {[G 0] + [sC -wC]}*[vd]
        % [iq]    [0 G]   [wC  sC]  [vq]
        if B == 0                           % G branch
            A_GC = []; B_GC = []; E_GC = [];
            C_GC = []; D_GC = inv([G,0;0,G]);
        else                                % GC or C branch
            % State equation
            % [d_vd] = 1/C*[-G  wC]*[vd] + 1/C*[1 0]*[id]
            % [d_vq]       [-wC -G] [vq]       [0 1] [iq]
            A_GC = 1/(B/w)*[-G,B;-B,-G];
            B_GC = 1/(B/w)*[1,0;0,1];
            E_GC = [1,0;0,1];
            % Output equation
            % [vd] = [1 0]*[vd] + [0 0]*[id]
            % [vq]   [0 1] [vq]   [0 0] [iq]
            C_GC = [1,0;0,1];
            D_GC = [0,0;0,0];
        end
        % Get the branch model
        Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);
        Y_GC= dss_SwitchInOut(Z_GC,2);
        Ybranch{n} = Y_GC;
    else                                % RL-GC blended branch, for special branch
        % RL branch
        if X == 0
            A_RL = []; B_RL = []; E_RL = []; C_RL = []; D_RL = inv([R,0;0,R]);
        else
            A_RL = 1/(X/w) * [-R,X;-X,-R];
            B_RL = 1/(X/w) * [1,0;0,1];
            E_RL = [1,0;0,1];
            C_RL = [1,0;0,1];
            D_RL = [0,0;0,0];
        end
        Y_RL = dss(A_RL,B_RL,C_RL,D_RL,E_RL);
        Z_RL = 1/Y_RL;

        % GC branch
        if B == 0
            A_GC = []; B_GC = []; E_GC = []; C_GC = []; D_GC = inv([G,0;0,G]);
        else
            A_GC = 1/(B/w)*[-G,B;-B,-G];
            B_GC = 1/(B/w)*[1,0;0,1];
            E_GC = [1,0;0,1];
            C_GC = [1,0;0,1];
            D_GC = [0,0;0,0];
        end
        Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);

        % Connect RL GC
        Z_RL_GC = Z_RL + Z_GC;
        Y_RL_GC = 1/Z_RL_GC;
        Ybranch{n} = Y_RL_GC;
    end

    % Get the state string of each branch
    if isempty(Ybranch{n}.A)
        Ybranch_StateStr{n}{1,1} = {''};
    else
        for i = 1:length(Ybranch{n}.A)
            Ybranch_StateStr{n}{1,i} = strcat('x_br',num2str(fb(n)),num2str(tb(n)),'_',num2str(i));
        end
    end
    
end

%%
% Initialize the state-space-form nodal addmitance matrix Ybus
Ass0 = []; Bss0 = []; Ess0 = [];
Css0 = []; Dss0 = [0,0;0,0];       % Defines a TITO static system 
Ybranch0 = dss(Ass0,Bss0,Css0,Dss0,Ess0);   % Descriptor state space model
for i = 1:n_bus
    for j = 1:n_bus
        YbusCell{i,j} = Ybranch0;
        YbusCell_StateStr{i,j} = {};
    end
end

% Calculate the element of the state-space-form nodal admittance matrix
for k=1:n_branch
    if fb(k) ~= tb(k) 
        % Off diagonal
        YbusCell{fb(k),tb(k)} = dss_Sum(YbusCell{fb(k),tb(k)},-Ybranch{k});
        YbusCell{tb(k),fb(k)} = dss_Sum(YbusCell{tb(k),fb(k)},-Ybranch{k});
        % Get state string
        YbusCell_StateStr{fb(k),tb(k)} = [YbusCell_StateStr{fb(k),tb(k)},Ybranch_StateStr{k}];
        YbusCell_StateStr{tb(k),fb(k)} = [YbusCell_StateStr{tb(k),fb(k)},Ybranch_StateStr{k}];
        
        % Diagonal
        YbusCell{fb(k),fb(k)} = dss_Sum(YbusCell{fb(k),fb(k)},Ybranch{k});
        YbusCell{tb(k),tb(k)} = dss_Sum(YbusCell{tb(k),tb(k)},Ybranch{k});
        % Get state string
     	YbusCell_StateStr{fb(k),fb(k)} = [YbusCell_StateStr{fb(k),fb(k)},Ybranch_StateStr{k}];
        YbusCell_StateStr{tb(k),tb(k)} = [YbusCell_StateStr{tb(k),tb(k)},Ybranch_StateStr{k}];
    else
        % Diagonal
        YbusCell{fb(k),tb(k)} = dss_Sum(YbusCell{fb(k),tb(k)},Ybranch{k});
        % Get state string
        YbusCell_StateStr{fb(k),tb(k)} = [YbusCell_StateStr{fb(k),tb(k)},Ybranch_StateStr{k}];
    end
end

%% Get the cell model
YbusCell;                           	% Cell-type state space form

%% Get the descriptor-state-space model
YbusDSS = dss_Arrange(YbusCell);      	% Whole state space form

%% Get the object model
% Create a new object
YbusObj = Class_Model_DSS;

% Load the model
YbusObj.LoadModel(YbusObj,YbusDSS);

% Get the string
for k = 1:n_bus
    InputStr{2*k-1}  = strcat('v_d',num2str(k));
    InputStr{2*k}    = strcat('v_q',num2str(k));
    OutputStr{2*k-1} = strcat('i_d',num2str(k));
    OutputStr{2*k}   = strcat('i_q',num2str(k));
end
StateStr = {};
for i = 1:n_bus
    for j = 1:n_bus
        % The sequence of connecting state strings is determined by "dss_Arange()"
        StateStr = [StateStr,YbusCell_StateStr{i,j}];
    end
end
YbusObj.WriteString(YbusObj,StateStr,InputStr,OutputStr);

% Check dimension mismatch
obj_CheckDim(YbusObj);

end

