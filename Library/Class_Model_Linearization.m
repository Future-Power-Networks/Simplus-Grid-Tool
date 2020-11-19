% This class defines the linearization algorithm of a general state space
% model:

% Author(s): Yitong Li, Yunjie Gu

%% Notes
% Format of original model
% dx/dt = f(x,u)
% y     = g(x,u)
% Format of linearized model
% dx/dt = A*x + B*u
% y     = C*x + D*u

%% Class
classdef Class_Model_Linearization < Class_Model_Base
    
    methods(Static)
        % Linearize state and output equations to get the linearized state
        % space matrices at a given steady-state operating point
        function Linearization(obj,x_e,u_e)

            % Calculate equilibrium of dx_e and y_e
            dx_e = obj.StateSpaceEqu(obj, x_e, u_e, 1);
            y_e  = obj.StateSpaceEqu(obj, x_e, u_e, 2);

            % Calculate length
            lx = length(x_e);
            lu = length(u_e);
            ly = length(y_e);

            % Initialize A,B,C,D
            A = zeros(lx,lx);
            B = zeros(lx,lu);
            C = zeros(ly,lx);
            D = zeros(ly,lu);

            % Get the perturb size
            perturb_factor = 1e-6;

            % Perturb x to calculate Ass and Css
            for i = 1:lx
                x_p = x_e;                      % Reset xp
                perturb_x = perturb_factor * abs(1+abs(x_e(i))); 	% Perturb size
                x_p(i) = x_e(i) + perturb_x;                        % Positive perturb on the ith element of xp
                dx_p = obj.StateSpaceEqu(obj, x_p, u_e, 1);
                y_p  = obj.StateSpaceEqu(obj, x_p, u_e, 2);
                A(:,i) = (dx_p - dx_e)/(x_p(i) - x_e(i));
                C(:,i) = (y_p - y_e)/(x_p(i) - x_e(i));
            end

            % Perturb u to calculate Bss and Dss
            for i = 1:lu
                up = u_e;                       % Reset up
                perturb_u = perturb_factor * abs(1+abs(u_e(i)));    % Perturb size
                up(i) = u_e(i) + perturb_u;                         % Positve perturb on the ith element of up
                dx_p = obj.StateSpaceEqu(obj, x_e, up, 1);
                y_p  = obj.StateSpaceEqu(obj, x_e, up, 2);
                B(:,i) = (dx_p - dx_e)/(up(i) - u_e(i));
                D(:,i) = (y_p - y_e)/(up(i) - u_e(i));
            end

            obj.MatrixSS    = {A,B,C,D};
            obj.A = A; obj.B = B; obj.C = C; obj.D = D;
        end
    end

end     % End class definition