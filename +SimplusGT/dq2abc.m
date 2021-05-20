% This function transforms a time-domain signal from the synchronous frame
% (dq0 frame) to natural frame (abc frame).

% Author(s): Yitong Li

% Notes:
%
% theta is the angle that dq axis leads abc axis.
%
% q axis leads d axis by 90 degrees.

function u_abc = dq2abc(u_dq,theta,varargin)

    % Default settings
    Flag_PowerVariant = false; % Default, power invariant
    Flag_DirectMatrixCalculation = true; % Default, direct calculation
    
    % Check advanced settings
    for n = 1:length(varargin)
        if (strcmpi(varargin{n},'PowerVariant'))
            Flag_PowerInvariant = varargin{n+1};
        elseif (strcmpi(varargin{n},'DirectMatrixCalculation'))
            Flag_DirectMatrixCalculation = varargin{n+1};
        end
    end
    
    % Check the form of the input signal
 	[r,c] = size(u_dq);
    if (r==2) && (c==1)
        Flag_ZeroSequence = 0;
    elseif (r==3) && (c==1)
        Flag_ZeroSequence = 1;
    else
        error('Error: Input signal is not a 2*1 or 3*1 vector.')
    end
    
    % Get the transformation matrix directly
    if Flag_DirectMatrixCalculation
        if Flag_PowerVariant
            T = [cos(theta),        -sin(theta),        1;
                 cos(theta-2/3*pi), -sin(theta-2/3*pi), 1;
                 cos(theta-4/3*pi), -sin(theta-4/3*pi), 1];
        else
            T = sqrt(2/3)*[cos(theta),        -sin(theta),        1/sqrt(2);
                           cos(theta-2/3*pi), -sin(theta-2/3*pi), 1/sqrt(2);
                           cos(theta-4/3*pi), -sin(theta-4/3*pi), 1/sqrt(2)];
        end
        
        % Check zero sequence
        if Flag_ZeroSequence
            T_dq2abc = T;
        else
            T_dq2abc = T(:,1:2);
        end
        
        % Frame transformation
        u_abc = T_dq2abc * u_dq;
        
  	% Call "abc2alphabeta" and "alphabeta2dq" to acheive two-step
    % transformation.
    else
        % Frame transformation
        u_alphabeta = SimplusGT.dq2alphabeta(u_dq,theta);
        u_abc = SimplusGT.alphabeta2abc(u_alphabeta,'PowerVariant',Flag_PowerVariant,'ZeroSequence',Flag_ZeroSequence);
    end

end