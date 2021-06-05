% This function calculate the descriptor (implicit) state space model of
% the dynamic nodal admittance matrix for ac networks.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% This function also considers the load impedance RL and XL.
%
% Self-branch: branches with from = to, connecting the bus and ground
% Mutual-branch: branches with from ~= to, connecting two difference buses
%
% The obtained nodal admittance matrix is in Laplace s domain, i.e., a
% dynamic model,is different from the steady nodal admittance matrix used
% in power flow analysis and only holds a similar form.
%
% This function can deal with both AC and DC grids. The default is AC. For
% AC grids, a two-input-two-output (TITO) transfer function or state space
% in dq frame for each branch; For DC grids, a single-input-single-output
% (SISO) transfer function or state space for each branch.
%
% For the format of the ListLine, please refer to comments in the
% "UserData".

%% Function

function [YbusObj,YbusDSS,YbusCell] = YbusCalcDss(ListBus,ListLine,w) 

%%
% Load the data
FB    = ListLine(:,1);              % From bus number...
TB    = ListLine(:,2);              % To bus number...
Rlist = ListLine(:,3);              % Resistance,  R...
Xlist = ListLine(:,4);              % Inductance,  wL...
Blist = ListLine(:,5);              % Capacitance, wC...
Glist = ListLine(:,6);              % Conductance, G...
Tlist = ListLine(:,7);              % Turns ratio, T

N_Bus = max(max(FB),max(TB));           % Number of buses
N_Branch = length(FB);                  % Number of branches, including self branches
      
% Inductive load effect
XLlist = ListLine(:,8);
AreaTypeLine = ListLine(:,9);
AreaType = ListBus(:,12);

%%
% Calculate the state space model of each branch
for n = 1:N_Branch              % Calculate the branch paramter one by one

    R = Rlist(n);
    X = Xlist(n);
    B = Blist(n);
    G = Glist(n);
    T = Tlist(n);
    XL = XLlist(n);

    if AreaTypeLine(n) == 1
    % ===========================
    % For AC grid
    % ============================
        if ( isinf(R) || isinf(X) || ( (G==0) && (B==0) && isinf(XL) ) )     % open circuit
            A_op = []; B_op = []; E_op = [];
            C_op = []; D_op = [0,0;0,0];
            Ybranch{n} = dss(A_op,B_op,C_op,D_op,E_op);
        elseif ( (R==0) && (X==0) && ( isinf(G) || isinf(B) || (XL==0) ) )
            error(['Error: short circuit, branch from ' num2str(FB(n)) ' to ' num2str(TB(n))]);
        else
            % ### Mutual branch
            if FB(n)~=TB(n)
                if ~( isinf(G) || isinf(B) )
                    error(['Error: Mutual branch ' num2str(FB(n)) num2str(TB(n)) ' contains G or B']);
                end
                if X == 0                       % R branch
                    A_RL = []; B_RL = []; E_RL = [];
                    C_RL = []; D_RL = inv([R,0;0,R]);
                else                            % RL or L branch
                    % KVL equations
                    % [vd] = {[R  0] + [sL -wL]}*[id]
                    % [vq]    [0  R]   [wL  sL]  [iq]
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
            % ### Self branch
            else
                % Error check
                if ~( (R==0) && (X==0) )         % GC branch, normally for self branch
                    error(['Error: Self branch ' num2str(FB(n)) num2str(TB(n)) ' contains R or X.']);
                end
                if B == 0                           % G branch
                    A_GC = []; B_GC = []; E_GC = [];
                    C_GC = []; D_GC = inv([G,0;0,G]);
                else                                % GC or C branch
                    % KCL equations
                    % [id] = {[G 0] + [sC -wC]}*[vd]
                    % [iq]    [0 G]   [wC  sC]  [vq]
                    % State equation
                    % [d_vd]/dt = 1/C*[-G  wC]*[vd] + 1/C*[1 0]*[id]
                    % [d_vq]          [-wC -G] [vq]       [0 1] [iq]
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
                Y_GC= SimplusGT.DssSwitchInOut(Z_GC,2);
                Ybranch{n} = Y_GC;
                % For self branch, connect load to it
                if isinf(XL)
                    Ybranch{n} = Ybranch{n};
                elseif (XL)==0
                    error(['Error: The inductive load is short-circuit. Please check QLi settings.']);
                else
                    % KVL equation for XL
                    % [vd] = {[sL -wL]}*[id]
                    % [vq]    [wL  sL]  [iq]
                    % => State equation
                    % [d_id]/dt = -1/L*[0 -wL]*[id] + 1/L*[1 0]*[vd] 
                    % [d_iq]           [wL  0] [iq]       [0 1] [vq]
                    A_XL = -1/(XL/w)*[0,-XL;XL,0];
                    B_XL = 1/(XL/w)*[1,0;0,1];
                    E_XL = [1,0;0,1];
                    % => Output equation
                    % [id] = [1 0]*[id] + [0 0]*[vd]
                    % [iq]   [0 1] [iq]   [0 0] [vq]
                    C_XL = [1,0;0,1];
                    D_XL = [0,0;0,0];

                    Y_XL = dss(A_XL,B_XL,C_XL,D_XL,E_XL);
                    Ybranch{n} = SimplusGT.DssSum(Ybranch{n},Y_XL);
                end
            end
        end
        
    elseif AreaTypeLine(n) == 2
  	% ===========================
    % For DC grid
    % ============================
        if ( isinf(R) || isinf(X) || ( (G==0) && (B==0) && isinf(XL) ) )     % open circuit
            A_op = []; B_op = []; E_op = [];
            C_op = []; D_op = [0];
            Ybranch{n} = dss(A_op,B_op,C_op,D_op,E_op);
        elseif ( (R==0) && (X==0) && ( isinf(G) || isinf(B) || (XL==0) ) )
            error(['Error: short circuit, branch from ' num2str(FB(n)) ' to ' num2str(TB(n))]);
        else
            % ### Mutual branch
            if FB(n)~=TB(n)
                if ~( isinf(G) || isinf(B) )
                    error(['Error: Mutual branch ' num2str(FB(n)) num2str(TB(n)) ' contains G or B']);
                end
                if X == 0                       % R branch
                    A_RL = []; B_RL = []; E_RL = [];
                    C_RL = []; D_RL = inv([R]);
                else                            % RL or L branch
                    % KVL equations
                    % v = (R + sL)*i
                    % State equations
                    % [d_i]/dt = 1/L*(-R)*[i] + 1/L*[v]
                    A_RL = 1/(X/w)*(-R);
                    B_RL = 1/(X/w);
                    E_RL = 1;
                    % Output equations
                    % [i] = 1*[i] + 0*[v]
                    C_RL = 1;
                    D_RL = 0;
                end
                % Get the branch model
                Y_RL = dss(A_RL,B_RL,C_RL,D_RL,E_RL); 
                Ybranch{n} = Y_RL;

            % ### Self branch
            else
                % Error check
                if ~( (R==0) && (X==0) )         % GC branch, normally for self branch
                    error(['Error: Self branch ' num2str(FB(n)) num2str(TB(n)) ' contains R or X.']);
                end
                if B == 0                           % G branch
                    A_GC = []; B_GC = []; E_GC = [];
                    C_GC = []; D_GC = inv([G]);
                else                                % GC or C branch
                    % KCL equations
                    % [i] = (G + sC)*[v]
                    % State equation
                    % [d_v]/dt = 1/C*[-G]*[v] + 1/C*[i]
                    A_GC = 1/(B/w)*(-G);
                    B_GC = 1/(B/w);
                    E_GC = 1;
                    % Output equation
                    % [v] = [1]*[v] + [0]*[i]
                    C_GC = 1;
                    D_GC = 0;
                end
                % Get the branch model
                Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);
                Y_GC= SimplusGT.DssSwitchInOut(Z_GC,1);
                Ybranch{n} = Y_GC;

                % For self branch, connect load to it
                if isinf(XL)
                    Ybranch{n} = Ybranch{n};
                elseif (XL)==0
                    error(['Error: The inductive load is short-circuit. Please check QLi settings.']);
                else
                    % KVL equation for XL
                    % [v] = [sL]*[i]
                    % => State equation
                    % d_i/dt = 0*[i] + 1/L*[v] 
                    A_XL = 0;
                    B_XL = 1/(XL/w);
                    E_XL = 1;
                    % => Output equation
                    % [i] = [1]*[i] + [0]*[v]
                    C_XL = 1;
                    D_XL = 0;

                    Y_XL = dss(A_XL,B_XL,C_XL,D_XL,E_XL);
                    Ybranch{n} = SimplusGT.DssSum(Ybranch{n},Y_XL);
                end
            end
        end
    
    end

    % Get the state string of each branch
    if isempty(Ybranch{n}.A)
        Ybranch_StateStr{n}{1,1} = {};
    else
        for i = 1:length(Ybranch{n}.A)
            Ybranch_StateStr{n}{1,i} = strcat('x_br',num2str(FB(n)),'-',num2str(TB(n)),'_',num2str(i));
        end
    end
    
end

%%
% Initialize the state-space-form nodal addmitance matrix Ybus
Ass0 = []; Bss0 = []; Ess0 = []; Css0 = [];
Dss0_ac2ac = [0,0;
              0,0];       % Defines a TITO static system for dq frame ac system.
Dss0_dc2dc = 0;               % Defines a SISO static system for dc system.
Dss0_ac2dc = [0,0];
Dss0_dc2ac = [0;
              0];

Ybr0_ac2ac = dss(Ass0,Bss0,Css0,Dss0_ac2ac,Ess0);
Ybr0_dc2dc = dss(Ass0,Bss0,Css0,Dss0_dc2dc,Ess0);
Ybr0_ac2dc = dss(Ass0,Bss0,Css0,Dss0_ac2dc,Ess0);
Ybr0_dc2ac = dss(Ass0,Bss0,Css0,Dss0_dc2ac,Ess0);

for i = 1:N_Bus
    for j = 1:N_Bus
        if (AreaType(i)==1) && (AreaType(j)==1)
            YbusCell{i,j} = Ybr0_ac2ac;
        elseif (AreaType(i)==2) && (AreaType(j)==2)
            YbusCell{i,j} = Ybr0_dc2dc;
      	elseif (AreaType(i)==1) && (AreaType(j)==2)
            YbusCell{i,j} = Ybr0_dc2ac;
        elseif (AreaType(i)==2) && (AreaType(j)==1)
            YbusCell{i,j} = Ybr0_ac2dc;
        else
            error(['Error']);
        end
        YbusCell_StateStr{i,j} = {};
    end
end

% Calculate the element of the state-space-form nodal admittance matrix
for k=1:N_Branch
    if FB(k) ~= TB(k) 
        % Off diagonal
        YbusCell{FB(k),TB(k)} = SimplusGT.DssSum(YbusCell{FB(k),TB(k)},-Ybranch{k}/Tlist(k));
        YbusCell{TB(k),FB(k)} = YbusCell{FB(k),TB(k)};
        % Get state string
        YbusCell_StateStr{FB(k),TB(k)} = [YbusCell_StateStr{FB(k),TB(k)},Ybranch_StateStr{k}];
        YbusCell_StateStr{TB(k),FB(k)} = YbusCell_StateStr{FB(k),TB(k)};
        
        % Diagonal
        YbusCell{FB(k),FB(k)} = SimplusGT.DssSum(YbusCell{FB(k),FB(k)},Ybranch{k}/Tlist(k)^2);
        YbusCell{TB(k),TB(k)} = SimplusGT.DssSum(YbusCell{TB(k),TB(k)},Ybranch{k});
        % Get state string
     	YbusCell_StateStr{FB(k),FB(k)} = [YbusCell_StateStr{FB(k),FB(k)},Ybranch_StateStr{k}];
        YbusCell_StateStr{TB(k),TB(k)} = [YbusCell_StateStr{TB(k),TB(k)},Ybranch_StateStr{k}];
    else
        % Diagonal
        YbusCell{FB(k),TB(k)} = SimplusGT.DssSum(YbusCell{FB(k),TB(k)},Ybranch{k});
        % Get state string
        YbusCell_StateStr{FB(k),TB(k)} = [YbusCell_StateStr{FB(k),TB(k)},Ybranch_StateStr{k}];
    end
end

%% Get the cell model
YbusCell;                           	% Cell-type state space form

%% Get the descriptor-state-space model
YbusDSS = SimplusGT.DssArrange(YbusCell);      	% Whole state space form

%% Get the object model
% Create a new object
YbusObj = SimplusGT.Class.ModelBase;

% Load the model
YbusObj.SetDSS(YbusObj,YbusDSS);

% Get the string
InOutCount = 1;
for k = 1:N_Bus
    if AreaType(k) == 1
        InputStr{InOutCount}    = strcat('v_d',num2str(k));
        InputStr{InOutCount+1}	= strcat('v_q',num2str(k));
        OutputStr{InOutCount}   = strcat('i_d',num2str(k));
        OutputStr{InOutCount+1}	= strcat('i_q',num2str(k));
        InOutCount = InOutCount + 2;
    elseif AreaType(k) == 2
        InputStr{InOutCount}    = strcat('v',num2str(k));
        OutputStr{InOutCount}	= strcat('i',num2str(k));
        InOutCount = InOutCount + 1;
    end
end
StateStr = {};
for i = 1:N_Bus
    for j = 1:N_Bus
        % The sequence of connecting state strings is determined by "dss_Arange()"
        StateStr = [StateStr,YbusCell_StateStr{i,j}];
    end
end
YbusObj.SetString(YbusObj,StateStr,InputStr,OutputStr);

% Check dimension mismatch
SimplusGT.ObjCheckDim(YbusObj);

end

