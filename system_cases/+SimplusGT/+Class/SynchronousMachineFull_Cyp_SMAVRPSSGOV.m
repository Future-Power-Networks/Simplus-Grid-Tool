% This class defines the model of 8-th order synchronous machine plus
% IEEE-T1 Exciter

% Author(s): Yue Zhu

%% Class

classdef SynchronousMachineFull_Cyp_SMAVRPSSGOV < SimplusGT.Class.ModelAdvance
    
    properties(Access = protected)
        psi_f;
        L;
        SE;
    end
    
    methods
        % constructor
        function obj = SynchronousMachineFull_Cyp_SMAVRPSSGOV(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        %% signal list
        
        function [State,Input,Output] = SignalList(obj)
                 % Vx1 and Vx2 are two middle states for PID.
          State = {'id','iq','w','theta','Ed1','Eq1','Psi1d','Psi2q',...
              'Efd','V1','Vr','Vf',...
              'Vssw','Vss1','Vss',...
              'Ptg1','Ptg2','Ptg3','T_m_in'};
          Input	 = {'v_d','v_q','Pmref','Vref'};
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
            Sbase_SM = obj.Para(1);
            X=obj.Para(2)/Sbase_SM;
            R=obj.Para(3)/Sbase_SM;
            Xd=obj.Para(4)/Sbase_SM; %synchronous reactance in d axis
            Xd1=obj.Para(5)/Sbase_SM; %transient reactance
            Xd2=obj.Para(6)/Sbase_SM; %subtransient reactance
            Td1=obj.Para(7); %d-axis open circuit transient time constant
            Td2=obj.Para(8); %d-axis open circuit sub-transient time constant
            Xq=obj.Para(9)/Sbase_SM;
            Xq1=obj.Para(10)/Sbase_SM;
            Xq2=obj.Para(11)/Sbase_SM;
            Tq1=obj.Para(12);
            Tq2=obj.Para(13);
            H=obj.Para(14)*Sbase_SM;
            D=obj.Para(15)*Sbase_SM;
            Dpu=obj.Para(16)*Sbase_SM;
            S10=obj.Para(17);
            S12=obj.Para(18);          
            ws=obj.Para(48);
            %AVR
            Tr=obj.Para(19);
            Ka=obj.Para(20);
            Ta=obj.Para(21);
            Vrmax=obj.Para(22);
            Vrmin=obj.Para(23);
            Ke=obj.Para(24);
            Te=obj.Para(25);
            Kf=obj.Para(26);
            Tf=obj.Para(27);
            E1=obj.Para(28);
            SEE1=obj.Para(29);
            E2=obj.Para(30);
            SEE2=obj.Para(31);
            %PSS
            %PSS
            T1=obj.Para(32);
            T2=obj.Para(33);
            T3=obj.Para(34);
            T4=obj.Para(35);
            Tw=obj.Para(36);
            Kpss=obj.Para(37);
            Vpssmin=obj.Para(38);
            Vpssmax=obj.Para(39);
            %Governor
            Rgov=obj.Para(40)/Sbase_SM;
            T1gov=obj.Para(41);
            T2gov=obj.Para(42);
            T3gov=obj.Para(43);
            T4gov=obj.Para(44);
            T5gov=obj.Para(45);
            Fgov=obj.Para(46);
            Pmax_gov = obj.Para(47)*Sbase_SM;

            % Calculate Equilibriums
            i_DQ = (conj(P+1j*Q)/V);
            Eint = V- i_DQ*(R+1j*Xq); %internal voltage
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
            Psi_q = R*i_d-v_d;
            Psi_d = -R*i_q+v_q;            
            T_e = Psi_d*i_q-Psi_q*i_d;%Psi_q*i_d-Psi_d*i_q;
            T_m_in = T_e - Dpu*w/ws;
            
            %PSS - conventional
            Vssw=0;
            Vss=0;
            Vss1=0;           
            % exciter IEEET1
            Vt=abs(v_d+1j*v_q);
            
            Bex=log(SEE1/SEE2)/(E1-E2);
            Aex=SEE1*exp(-Bex*E1);
                      
            Vss=0;
            Vf=0;          
            V1=Vt;
            Vr=Ke*Efd+Efd*Aex*exp(Bex*Efd);
            if Vr>=Vrmax
                Vr=Vrmax;
            elseif Vr<=Vrmin
                Vr=Vrmin;
            end      
            Vref = Vr/Ka+Vf+V1-Vss;
            
            % BPA-GG governor
            
            Ptg1 = 0;
            Ptg3 = T_m_in+2*Sbase_SM; % see the figure in power factory.           
            % if Ptg3<=-Pmax_gov
            %     Ptg3 = -Pmax_gov;
            % elseif Ptg3>=0
            %     Ptg3 = 0;
            % end            
            Ptg2 = Ptg3;
            Pmref = Ptg2;
                        
            % angle rotate
            xi = xi+sigma-pi/2;
            theta = xi;
                           
            x_e = [i_d; i_q; w; theta; Ed1; Eq1; Psi1d; Psi2q; ...
                Efd; V1; Vr; Vf;...
                Vssw; Vss1; Vss;...
                Ptg1; Ptg2; Ptg3; T_m_in];
            u_e = [v_d; v_q; Pmref; Vref];
           
            % Get equilibrium
            xi  = [xi];    %output of this function.        
         
        end
        
        %% State-space
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
                        % Get parameters
            % Get parameters
            % Synchronous machine
            Sbase_SM = obj.Para(1);
            X=obj.Para(2)/Sbase_SM;
            R=obj.Para(3)/Sbase_SM;
            Xd=obj.Para(4)/Sbase_SM; %synchronous reactance in d axis
            Xd1=obj.Para(5)/Sbase_SM; %transient reactance
            Xd2=obj.Para(6)/Sbase_SM; %subtransient reactance
            Td1=obj.Para(7); %d-axis open circuit transient time constant
            Td2=obj.Para(8); %d-axis open circuit sub-transient time constant
            Xq=obj.Para(9)/Sbase_SM;
            Xq1=obj.Para(10)/Sbase_SM;
            Xq2=obj.Para(11)/Sbase_SM;
            Tq1=obj.Para(12);
            Tq2=obj.Para(13);
            H=obj.Para(14)*Sbase_SM;
            D=obj.Para(15)*Sbase_SM;
            Dpu=obj.Para(16)*Sbase_SM;
            S10=obj.Para(17);
            S12=obj.Para(18);          
            ws=obj.Para(48);
            %AVR
            Tr=obj.Para(19);
            Ka=obj.Para(20);
            Ta=obj.Para(21);
            Vrmax=obj.Para(22);
            Vrmin=obj.Para(23);
            Ke=obj.Para(24);
            Te=obj.Para(25);
            Kf=obj.Para(26);
            Tf=obj.Para(27);
            E1=obj.Para(28);
            SEE1=obj.Para(29);
            E2=obj.Para(30);
            SEE2=obj.Para(31);
            %PSS
            %PSS
            T1=obj.Para(32);
            T2=obj.Para(33);
            T3=obj.Para(34);
            T4=obj.Para(35);
            Tw=obj.Para(36);
            Kpss=obj.Para(37);
            Vpssmin=obj.Para(38);
            Vpssmax=obj.Para(39);
            %Governor
            Rgov=obj.Para(40)/Sbase_SM;
            T1gov=obj.Para(41);
            T2gov=obj.Para(42);
            T3gov=obj.Para(43);
            T4gov=obj.Para(44);
            T5gov=obj.Para(45);
            Fgov=obj.Para(46);
            Pmax_gov = obj.Para(47)*Sbase_SM;
            
            i_d = x(1);
            i_q = x(2);
            w = x(3);
            theta = x(4);
            Ed1 = x(5);
            Eq1 = x(6);
            Psi1d = x(7);
            Psi2q = x(8);  
            Efd = x(9);
            V1 = x(10);
            Vr = x(11);
            Vf = x(12);
            Vssw=x(13);
            Vss1=x(14);
            Vss=x(15);
            Ptg1=x(16);
            Ptg2=x(17);
            Ptg3=x(18);
            T_m_in=x(19);
                
            v_d  = u(1);
            v_q  = u(2);
            Pmref  = u(3);
            Vref  = u(4);
            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)

            if CallFlag == 1        
            % ### Call state equation: dx/dt = f(x,u)
                % Auxiliary equations
                Psi_d = Xd2*i_d+(Xd2-X)/(Xd1-X)*Eq1 + (Xd1-Xd2)/(Xd1-X)*Psi1d; %book 
                Psi_q = Xq2*i_q-(Xq2-X)/(Xq1-X)*Ed1 + (Xq1-Xq2)/(Xq1-X)*Psi2q;
                T_e = Psi_d*i_q-Psi_q*i_d;%T_e = (Psi_q*i_d-Psi_d*i_q);
                % SM
                dtheta = w;
                dw     = ws*(T_e - T_m_in - D*(w-ws)/ws - Dpu*w/ws)/(2*H); % original
                dEd1 = (-Ed1+(Xq-Xq1)*(-i_q+(Xq1-Xq2)/(Xq1-X)^2*(-Psi2q+(Xq1-X)*i_q-Ed1)))/Tq1;
                dEq1 = (Efd-Eq1+(Xd-Xd1)*(i_d+(Xd1-Xd2)/(Xd1-X)^2*(Psi1d-(Xd1-X)*i_d-Eq1)))/Td1;                 
                dPsi1d = (-Psi1d+Eq1+(Xd1-X)*i_d)/Td2;
                dPsi2q = (-Psi2q-Ed1+(Xq1-X)*i_q)/Tq2;
                di_d = (ws*(v_d-R*i_d+w/ws*Psi_q)-(Xd2-X)/(Xd1-X)*dEq1-(Xd1-Xd2)/(Xd1-X)*dPsi1d)/Xd2;
                di_q = (ws*(v_q-R*i_q-w/ws*Psi_d)+(Xq2-X)/(Xq1-X)*dEd1-(Xq1-Xq2)/(Xq1-X)*dPsi2q)/Xq2;
                % PSS  state equations
                if Vss>=Vpssmax
                    Vss=Vpssmax;
                elseif Vss<=Vpssmin
                    Vss=Vpssmin;
                end
                dVssw = (Kpss*Tw*dw-Vssw)/Tw;
                dVss1 = (Vssw+T1*dVssw-Vss1)/T2;
                dVss  = (Vss1+T3*dVss1-Vss)/T4;                         
                % Exciter state equations
                if Vr>=Vrmax
                    Vr=Vrmax;
                elseif Vr<=Vrmin
                    Vr=Vrmin;
                end
                Bex=log(SEE1/SEE2)/(E1-E2);
                Aex=SEE1*exp(-Bex*E1);
                Vt=abs(v_d+1j*v_q);
                dEfd = (Vr-(Ke*Efd+Efd*Aex*exp(Bex*Efd)))/Te;
                if Tr~=0
                    dV1=(Vt-V1)/Tr;
                else
                    dV1 = 0;
                    V1 = Vt;
                end
                dVr = ((Vss-V1+Vref-Vf)*Ka-Vr)/Ta;
                dVf = (Kf*dEfd-Vf)/Tf;
                
                % Governor state equations
                K=1/Rgov;
                Delta_w=-(ws-w)/ws;
                if T1gov~=0
                    dPtg1 = (K*Delta_w-K*T2gov*dw-Ptg1)/T1gov;
                else
                    dPtg1 = 0;%-K*dw;
                    Ptg1 = K*Delta_w;
                end
                %dPtg1 = ( (Delta_w-T2gov*dw)/Rgov-Ptg1 ) / T1gov;
                Px=Ptg1+Pmref;
                % if Px<=-Pmax_gov
                %     Px=-Pmax_gov;
                % elseif Px>=0
                %     Px=0;
                % end
                dPtg2 = (Px-Ptg2)/T3gov;
                dPtg3 = (Ptg2-Ptg3)/T4gov;
                %if T5gov~=0
                if T5gov == 0
                    dT_m_in = 0;
                    %T_m_in = Ptg3;
                    %T5gov = 1e-7;
                %end
                else
                    dT_m_in  = (Ptg3+Fgov*T5gov*dPtg3-T_m_in-2*Sbase_SM)/(T5gov);
                end               
                % end of state calculation
                f_xu = [di_d; di_q; dw; dtheta; dEd1; dEq1; dPsi1d; dPsi2q;...
                    dEfd; dV1; dVr; dVf;...
                    dVssw; dVss1; dVss;...
                    dPtg1; dPtg2; dPtg3; dT_m_in];
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