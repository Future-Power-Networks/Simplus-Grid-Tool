% This function transforms a time-domain signal from synchronous frame (dq0
% frame) to stationary frame (alpha/beta/0 frame).

% Author(s): Yitong Li

% Notes:
%
% theta is the angle that dq axis leads alpha/beta axis.
%
% q axis leads d axis by 90 degrees.

function u_alphabeta = dq2alphabeta(u_dq,theta)

    % Check the form of input signal
    [r,c] = size(u_dq);

    if (r==2) && (c==1)
        Flag_ZeroSequence = 0;
    elseif (r==3) && (c==1)
        Flag_ZeroSequence = 1;
    else
        error('Error: Input signal is not a 2*1 or 3*1 vector.')
    end
    
    % Get the inverse Park transformation matrix
 	T_inv_Park = [cos(theta), -sin(theta),0;
                  sin(theta),cos(theta),0;
                  0,          0,         1];
              
 	% Check zero sequence
    if Flag_ZeroSequence
        T_dq2alphabeta = T_inv_Park;
    else
        T_dq2alphabeta = T_inv_Park(1:2,1:2);
    end

    % Frame transformation
    u_alphabeta = T_dq2alphabeta * u_dq;

end