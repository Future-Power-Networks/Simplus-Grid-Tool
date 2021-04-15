% This class defines the model of 8-th order synchronous machine

% Author(s): Yue Zhu
% Reference Book: Power System Dynamic and Stability, P54

%% Class

classdef SynchronousMachineFull_SM < SimplexPS.Class.ModelAdvance
    
    properties(Access = protected)
        psi_f;
        ws;
        L;
    end
    
    methods
        % constructor
        function obj = SynchronousMachineFull_SM(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        %% signal list
        
        function [State,Input,Output] = SignalList(obj)
                State = {'id', 'iq', 'w', 'theta','Ed1','Eq1','Psi1d','Psi2q'};
                % theta: rotor angle | w: rotor velocity 
                % Ed1: transient emf due to flux linkage in d axis
                % Eq1: transient emf due to flux linkage in q axis
                % Psi1d: sub-transient emf due to d-axis damper coils 
                % Psi2q: sub-transient emf due to q-axis damper coils
                Input	 = {'v_d','v_q','T_m','Efd'};
                Output = {'i_d','i_q','w','theta'};
        end
        %% Equilibrium point
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);            
            % Get parameters
            % Synchronous machine
            X=obj.Para(1);
            R=obj.Para(2);
            Xd=obj.Para(3); %synchronous reactance in d axis
            Xd1=obj.Para(4); %transient reactance
            Xd2=obj.Para(5); %subtransient reactance
            Td1=obj.Para(6); %d-axis open circuit transient time constant
            Td2=obj.Para(7); %d-axis open circuit sub-transient time constant
            Xq=obj.Para(8);
            Xq1=obj.Para(9);
            Xq2=obj.Para(10);
            Tq1=obj.Para(11);
            Tq2=obj.Para(12);
            H=obj.Para(13);
            D=obj.Para(14);
            
            % Calculate Equilibriums
            i_DQ = (conj(P+1j*Q)/V);
            Eint = V - i_DQ*(R+1j*Xq); %internal voltage
            sigma = angle(Eint);           
            i_dq = i_DQ * exp(1j*(-sigma+pi/2));
            i_d = real(i_dq);
            i_q = imag(i_dq);
            v_dq = V * exp(1j*(-sigma+pi/2));
            v_d = real(v_dq);
            v_q = imag(v_dq);           
            Efd = abs(Eint)+(Xq-Xd)*i_d;
            Eq1 = Efd+(Xd-Xd1)*i_d;           
            Ed1 = -(Xq-Xq1)*i_q;
            Psi1d=Eq1+(Xd1-X)*i_d;
            Psi2q=-Ed1+(Xq1-X)*i_q;

            obj.ws = w; % record synchronous rotor speed.
            obj.L = X/w; % record inductor value

            Psi_q = R*i_d-v_d;
            Psi_d = -R*i_q+v_q;            
            T_e = Psi_d*i_q-Psi_q*i_d;%Psi_q*i_d-Psi_d*i_q;
            T_m = T_e;
            
            xi = xi+sigma-pi/2;
            theta = xi;
           
                % Get equilibrium
                x_e = [i_d; i_q; w; theta; Ed1; Eq1; Psi1d; Psi2q];
                u_e = [v_d; v_q; T_m; Efd];
            % Get equilibrium
            xi  = [xi];    %connection angle to the whole grid.

            
            
        end
        
        %% State-space
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
                        % Get parameters
            L=obj.L;
            ws=obj.ws;
            X=obj.Para(1);
            R=obj.Para(2);
            Xd=obj.Para(3); %synchronous reactance in d axis
            Xd1=obj.Para(4); %transient reactance
            Xd2=obj.Para(5); %subtransient reactance
            Td1=obj.Para(6); %d-axis open circuit transient time constant
            Td2=obj.Para(7); %d-axis open circuit sub-transient time constant
            Xq=obj.Para(8);
            Xq1=obj.Para(9);
            Xq2=obj.Para(10);
            Tq1=obj.Para(11);
            Tq2=obj.Para(12);
            H=obj.Para(13);
            D=obj.Para(14);

            % Get states
            i_d = x(1);
            i_q = x(2);
            w = x(3);
            theta = x(4);
            Ed1 = x(5);
            Eq1 = x(6);
            Psi1d = x(7);
            Psi2q = x(8);
            
            % Get input signals
            v_d  = u(1);
            v_q  = u(2);
            T_m  = u(3);
            Efd  = u(4);
                       
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1        
            % ### Call state equation: dx/dt = f(x,u)
                % Auxiliary equations
                Psi_d = Xd2*i_d+(Xd2-X)/(Xd1-X)*Eq1 + (Xd1-Xd2)/(Xd1-X)*Psi1d; %book 
                Psi_q = Xq2*i_q-(Xq2-X)/(Xq1-X)*Ed1 + (Xq1-Xq2)/(Xq1-X)*Psi2q;               
                ws=obj.ws;
                T_e = Psi_d*i_q-Psi_q*i_d;%T_e = (Psi_q*i_d-Psi_d*i_q);
                % State equations
                dtheta = w;
                dw     = ws*(T_e - T_m - (D/ws)*(w-ws))/(2*H);
                dEd1 = (-Ed1+(Xq-Xq1)*(-i_q+(Xq1-Xq2)/(Xq1-X)^2*(-Psi2q+(Xq1-X)*i_q-Ed1)))/Tq1;
                dEq1 = (Efd-Eq1+(Xd-Xd1)*(i_d+(Xd1-Xd2)/(Xd1-X)^2*(Psi1d-(Xd1-X)*i_d-Eq1)))/Td1;                 
                dPsi1d = (-Psi1d+Eq1+(Xd1-X)*i_d)/Td2;
                dPsi2q = (-Psi2q-Ed1+(Xq1-X)*i_q)/Tq2;
                di_d = (ws*(v_d-R*i_d+w/ws*Psi_q)-(Xd2-X)/(Xd1-X)*dEq1-(Xd1-Xd2)/(Xd1-X)*dPsi1d)/Xd2; 
                di_q = (ws*(v_q-R*i_q-w/ws*Psi_d)+(Xq2-X)/(Xq1-X)*dEd1-(Xq1-Xq2)/(Xq1-X)*dPsi2q)/Xq2;
                
                f_xu = [di_d; di_q; dw; dtheta; dEd1; dEq1; dPsi1d; dPsi2q];
         
                Output = f_xu;
            elseif CallFlag == 2    
            % ### Call output equation: y = g(x,u)
                %i_ex = 0;
                g_xu = [i_d; i_q; w; theta];
                Output = g_xu;
            end
        end
        
    end

end     % End class definition