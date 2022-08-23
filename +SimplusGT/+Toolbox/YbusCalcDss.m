% This function calculate the descriptor (implicit) state space model of
% the dynamic nodal admittance matrix for ac networks.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% This function also considers the passive LCG load.
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

NumBus = max(max(FB),max(TB));           % Number of buses
NumBr = length(FB);                  % Number of branches, including self branches
      
AreaTypeLine = ListLine(:,8);
AreaTypeBus = ListBus(:,12);

%%
% Calculate the state space model of each branch
for n = 1:NumBr              % Calculate the branch paramter one by one

    R = Rlist(n);
    X = Xlist(n);
    B = Blist(n);
    G = Glist(n);
    T = Tlist(n);

    if AreaTypeLine(n) == 1
    % ===========================
    % For AC grid
    % ============================
        if (isinf(R) || isinf(X)) && (G==0) && (B==0)      % open circuit
            A_op = []; B_op = []; E_op = [];
            C_op = []; D_op = [0,0;0,0];
            YbrDss{n} = dss(A_op,B_op,C_op,D_op,E_op);
            YbrStateStr{n} = {};
        elseif ( (R==0) && (X==0) ) || isinf(G) || isinf(B)
            error(['Error: short circuit, ac branch from ' num2str(FB(n)) ' to ' num2str(TB(n))]);
        else
            % ### Mutual branch
            if FB(n)~=TB(n)
                if ~( G==0 && B==0 )
                    error(['Error: Ac mutual branch ' num2str(FB(n)) num2str(TB(n)) ' contains C and/or G']);
                end
                if X == 0                       % R branch
                    A_RL = []; B_RL = []; E_RL = [];
                    C_RL = []; D_RL = inv([R,0;0,R]);
                    YbrStateStr{n} = {};
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
                    YbrStateStr{n}(1) = {['id',num2str(FB(n)),'-',num2str(TB(n))]};
                    YbrStateStr{n}(2) = {['iq',num2str(FB(n)),'-',num2str(TB(n))]};
                end
                % Get the branch model
                Y_RL = dss(A_RL,B_RL,C_RL,D_RL,E_RL); 
                YbrDss{n} = Y_RL;
            % ### Self branch
            else
                % Error check
                if (R~=0) && (~isinf(R))            % LCG branch, normally for self branch
                    error(['Error: Ac self branch ' num2str(FB(n)) num2str(TB(n)) ' contains R.']);
                end
                if B == 0                           % G branch
                    A_GC = []; B_GC = []; E_GC = [];
                    C_GC = []; D_GC = inv([G,0;0,G]);
                    YbrStateStr{n} = {};
                else                                % CG or C branch
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
                    YbrStateStr{n}(1) = {['vd',num2str(FB(n)),'-',num2str(TB(n))]};
                    YbrStateStr{n}(2) = {['vq',num2str(FB(n)),'-',num2str(TB(n))]};
                end
                % Get the branch model
                Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);
                Y_GC= SimplusGT.DssSwitchInOut(Z_GC,2);
                YbrDss{n} = Y_GC;
                YbrStateStr{n} = [ YbrStateStr{n},...
                                   {['xi',num2str(FB(n)),'-',num2str(TB(n))]},...
                                   {['xi',num2str(FB(n)),'-',num2str(TB(n))]}];
                % For self branch, connect inductive load to it
                if isinf(X)
                    YbrDss{n} = YbrDss{n};
                elseif X==0
                    error(['Error: The inductive load is short-circuit. Please check QLi settings.']);
                else
                    % KVL equation for X
                    % [vd] = {[sL -wL]}*[id]
                    % [vq]    [wL  sL]  [iq]
                    % => State equation
                    % [d_id]/dt = -1/L*[0 -wL]*[id] + 1/L*[1 0]*[vd] 
                    % [d_iq]           [wL  0] [iq]       [0 1] [vq]
                    A_X = -1/(X/w)*[0,-X;X,0];
                    B_X = 1/(X/w)*[1,0;0,1];
                    E_X = [1,0;0,1];
                    % => Output equation
                    % [id] = [1 0]*[id] + [0 0]*[vd]
                    % [iq]   [0 1] [iq]   [0 0] [vq]
                    C_X = [1,0;0,1];
                    D_X = [0,0;0,0];

                    Y_X = dss(A_X,B_X,C_X,D_X,E_X);
                    YbrDss{n} = SimplusGT.DssSum(YbrDss{n},Y_X);
                    
                    YbrStateStr{n} = [ YbrStateStr{n},...
                                       {['id',num2str(FB(n)),'-',num2str(TB(n))]},...
                                       {['iq',num2str(FB(n)),'-',num2str(TB(n))]}];
                end
            end
        end
        
    elseif AreaTypeLine(n) == 2
  	% ===========================
    % For DC grid
    % ============================
        if ( isinf(R) || isinf(X) ) && (G==0) && (B==0)    % open circuit
            A_op = []; B_op = []; E_op = [];
            C_op = []; D_op = [0];
            YbrDss{n} = dss(A_op,B_op,C_op,D_op,E_op);
            YbrStateStr{n} = {};
        elseif ((R==0) && (X==0) ) || isinf(G) || isinf(B)
            error(['Error: short circuit, dc branch from ' num2str(FB(n)) ' to ' num2str(TB(n))]);
        else
            % ### Mutual branch
            if FB(n)~=TB(n)
                if ~( G==0 && B==0 )
                    error(['Error: Dc mutual branch ' num2str(FB(n)) '-' num2str(TB(n)) ' contains C and/or G']);
                end
                if X == 0                       % R branch
                    A_RL = []; B_RL = []; E_RL = [];
                    C_RL = []; D_RL = inv([R]);
                    YbrStateStr{n} = {};
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
                    YbrStateStr{n} = {['i',num2str(FB(n)),'-',num2str(TB(n))]};
                end
                % Get the branch model
                Y_RL = dss(A_RL,B_RL,C_RL,D_RL,E_RL); 
                YbrDss{n} = Y_RL;

            % ### Self branch
            else
                % Error check
                if ~( isinf(R) || isinf(X) )     	% CG branch, normally for self branch
                    error(['Error: Dc self branch ' num2str(FB(n)) num2str(TB(n)) ' contains R and/or L.']);
                end
                if B == 0                           % G branch
                    A_GC = []; B_GC = []; E_GC = [];
                    C_GC = []; D_GC = inv([G]);
                    YbrStateStr{n} = {};
                else                                % CG or C branch
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
                    YbrStateStr{n} = {['v',num2str(FB(n)),'-',num2str(TB(n))]};
                end
                % Get the branch model
                Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);
                Y_GC= SimplusGT.DssSwitchInOut(Z_GC,1);
                YbrStateStr{n} = [YbrStateStr{n},...
                                  {['xi',num2str(FB(n)),'-',num2str(TB(n))]}];
                YbrDss{n} = Y_GC;
            end
        end
    
    end
    
end

%%
% Initialize the state-space-form nodal addmitance matrix Ybus
Ass0 = []; Bss0 = []; Ess0 = []; Css0 = [];
Dss0_ac2ac = [0,0;
              0,0];         % Define a TITO static system for dq frame ac system.
Dss0_dc2dc = 0;            	% Define a SISO static system for dc system.
Dss0_ac2dc = [0,0];
Dss0_dc2ac = [0;
              0];

Ybr0_ac2ac = dss(Ass0,Bss0,Css0,Dss0_ac2ac,Ess0);
Ybr0_dc2dc = dss(Ass0,Bss0,Css0,Dss0_dc2dc,Ess0);
Ybr0_ac2dc = dss(Ass0,Bss0,Css0,Dss0_ac2dc,Ess0);
Ybr0_dc2ac = dss(Ass0,Bss0,Css0,Dss0_dc2ac,Ess0);

for i = 1:NumBus
    for j = 1:NumBus
        if (AreaTypeBus(i)==1) && (AreaTypeBus(j)==1)
            YbusCell{i,j} = Ybr0_ac2ac;
        elseif (AreaTypeBus(i)==2) && (AreaTypeBus(j)==2)
            YbusCell{i,j} = Ybr0_dc2dc;
      	elseif (AreaTypeBus(i)==1) && (AreaTypeBus(j)==2)
            YbusCell{i,j} = Ybr0_dc2ac;
        elseif (AreaTypeBus(i)==2) && (AreaTypeBus(j)==1)
            YbusCell{i,j} = Ybr0_ac2dc;
        else
            error(['Error']);
        end
        YbusCellStateStr{i,j} = {};
    end
end

% Calculate the element of the state-space-form nodal admittance matrix
for k=1:NumBr
    if FB(k) ~= TB(k) 
        % Off diagonal
        YbusCell{FB(k),TB(k)} = SimplusGT.DssSum(YbusCell{FB(k),TB(k)},-YbrDss{k}/Tlist(k));
        YbusCell{TB(k),FB(k)} = YbusCell{FB(k),TB(k)};
        % Get state string
        YbusCellStateStr{FB(k),TB(k)} = [YbusCellStateStr{FB(k),TB(k)},YbrStateStr{k}];
        YbusCellStateStr{TB(k),FB(k)} = YbusCellStateStr{FB(k),TB(k)};
        
        % Diagonal
        YbusCell{FB(k),FB(k)} = SimplusGT.DssSum(YbusCell{FB(k),FB(k)},YbrDss{k}/Tlist(k)^2);
        YbusCell{TB(k),TB(k)} = SimplusGT.DssSum(YbusCell{TB(k),TB(k)},YbrDss{k});
        % Get state string
     	YbusCellStateStr{FB(k),FB(k)} = [YbusCellStateStr{FB(k),FB(k)},YbrStateStr{k}];
        YbusCellStateStr{TB(k),TB(k)} = [YbusCellStateStr{TB(k),TB(k)},YbrStateStr{k}];
    else
        % Diagonal
        YbusCell{FB(k),TB(k)} = SimplusGT.DssSum(YbusCell{FB(k),TB(k)},YbrDss{k});
        % Get state string
        YbusCellStateStr{FB(k),TB(k)} = [YbusCellStateStr{FB(k),TB(k)},YbrStateStr{k}];
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
for k = 1:NumBus
    if AreaTypeBus(k) == 1
        InputStr{InOutCount}    = ['v_d',num2str(k)];
        InputStr{InOutCount+1}	= ['v_q',num2str(k)];
        OutputStr{InOutCount}   = ['i_d',num2str(k)];
        OutputStr{InOutCount+1}	= ['i_q',num2str(k)];
        InOutCount = InOutCount + 2;
    elseif AreaTypeBus(k) == 2
        InputStr{InOutCount}    = ['v',num2str(k)];
        OutputStr{InOutCount}	= ['i',num2str(k)];
        InOutCount = InOutCount + 1;
    end
end
StateStr = {};
for i = 1:NumBus
    for j = 1:NumBus
        % The order of connecting state strings is determined by "DssArange()"
        StateStr = [StateStr,YbusCellStateStr{i,j}];
    end
end
YbusObj.SetString(YbusObj,StateStr,InputStr,OutputStr);

% Check dimension mismatch
SimplusGT.ObjCheckDim(YbusObj);

end

