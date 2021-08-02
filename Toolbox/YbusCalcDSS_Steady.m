% This function calculate the steady state space model of nodal admittance
% matrix. The EMT dyanmics are ignored, for example, an inductor (s+jw)*L
% -> jw*L.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% This function also considers the load impedance RL and XL.
%
% Self-branch: branches with from = to, connecting the bus and ground
% Mutual-branch: branches with from ~= to, connecting two difference buses
%
% FrameFlag: 1-abc, 2-alpha/beta, 3-dq

%% Function

function [YbusObj,YbusDSS,YbusCell] = YbusCalcDSS(ListLine,w) 

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

%%
% Calculate the state space model of each branch
for n = 1:N_Branch              % Calculate the branch paramter one by one

    R = Rlist(n);
    X = Xlist(n);
    B = Blist(n);
    G = Glist(n);
    T = Tlist(n);
    XL = XLlist(n);

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
                % [vd] = {[R  0] + [0 -wL]}*[id]
                % [vq]    [0  R]   [wL  0]  [iq]
                A_RL = []; B_RL = []; E_RL = []; C_RL = [];
                D_RL = inv([R,-X;X,R]);
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
                % [id] = {[G 0] + [0  -wC]}*[vd]
                % [iq]    [0 G]   [wC   0]  [vq]
                A_GC = []; B_GC = []; E_GC = []; C_GC = [];
                D_GC = inv([G,-B;B,G]);
            end
            % Get the branch model
            Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);
            Y_GC= inv(Z_GC);
            Ybranch{n} = Y_GC;

            % For self branch, connect load to it
            if isinf(XL)
                Ybranch{n} = Ybranch{n};
            elseif (XL)==0
                error(['Error: The inductive load is short-circuit. Please check QLi settings.']);
            else
                % KVL equation for XL
                % [vd] = {[0  -wL]}*[id]
                % [vq]    [wL   0]  [iq]
                A_XL = []; B_XL = []; E_XL = []; C_XL = []; 
                D_XL = inv([0,-X;X,0]);

                Y_XL = dss(A_XL,B_XL,C_XL,D_XL,E_XL);
                Ybranch{n} = dss_Sum(Ybranch{n},Y_XL);
            end
        end
        
    end

    % Get the state string of each branch
    if isempty(Ybranch{n}.A)
        Ybranch_StateStr{n}{1,1} = {};
    else
        for i = 1:length(Ybranch{n}.A)
            Ybranch_StateStr{n}{1,i} = strcat('x_br',num2str(FB(n)),num2str(TB(n)),'_',num2str(i));
        end
    end
    
end

%%
% Initialize the state-space-form nodal addmitance matrix Ybus
Ass0 = []; Bss0 = []; Ess0 = [];
Css0 = []; Dss0 = [0,0;0,0];       % Defines a TITO static system 
Ybranch0 = dss(Ass0,Bss0,Css0,Dss0,Ess0);
for i = 1:N_Bus
    for j = 1:N_Bus
        YbusCell{i,j} = Ybranch0;
        YbusCell_StateStr{i,j} = {};
    end
end

% Calculate the element of the state-space-form nodal admittance matrix
for k=1:N_Branch
    if FB(k) ~= TB(k) 
        % Off diagonal
        YbusCell{FB(k),TB(k)} = dss_Sum(YbusCell{FB(k),TB(k)},-Ybranch{k}/Tlist(k));
        YbusCell{TB(k),FB(k)} = YbusCell{FB(k),TB(k)};
        % Get state string
        YbusCell_StateStr{FB(k),TB(k)} = [YbusCell_StateStr{FB(k),TB(k)},Ybranch_StateStr{k}];
        YbusCell_StateStr{TB(k),FB(k)} = YbusCell_StateStr{FB(k),TB(k)};
        
        % Diagonal
        YbusCell{FB(k),FB(k)} = dss_Sum(YbusCell{FB(k),FB(k)},Ybranch{k}/Tlist(k)^2);
        YbusCell{TB(k),TB(k)} = dss_Sum(YbusCell{TB(k),TB(k)},Ybranch{k});
        % Get state string
     	YbusCell_StateStr{FB(k),FB(k)} = [YbusCell_StateStr{FB(k),FB(k)},Ybranch_StateStr{k}];
        YbusCell_StateStr{TB(k),TB(k)} = [YbusCell_StateStr{TB(k),TB(k)},Ybranch_StateStr{k}];
    else
        % Diagonal
        YbusCell{FB(k),TB(k)} = dss_Sum(YbusCell{FB(k),TB(k)},Ybranch{k});
        % Get state string
        YbusCell_StateStr{FB(k),TB(k)} = [YbusCell_StateStr{FB(k),TB(k)},Ybranch_StateStr{k}];
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
for k = 1:N_Bus
    InputStr{2*k-1}  = strcat('v_d',num2str(k));
    InputStr{2*k}    = strcat('v_q',num2str(k));
    OutputStr{2*k-1} = strcat('i_d',num2str(k));
    OutputStr{2*k}   = strcat('i_q',num2str(k));
end
StateStr = {};
% for i = 1:N_Bus
%     for j = 1:N_Bus
%         % The sequence of connecting state strings is determined by "dss_Arange()"
%         StateStr = [StateStr,YbusCell_StateStr{i,j}];
%     end
% end

YbusObj.WriteString(YbusObj,StateStr,InputStr,OutputStr);

% Check dimension mismatch
obj_CheckDim(YbusObj);

end

