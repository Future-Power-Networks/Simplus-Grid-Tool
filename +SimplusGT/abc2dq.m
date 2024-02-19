% This function transforms a time-domain signal from the natural frame (abc
% frame) to synchronous frame (dq0 frame).

% Author(s): Yitong Li

% Notes:
%
% theta is the angle that dq axis leads abc axis.
% 
% q axis leads d axis by 90 degrees.

function u_dq = abc2dq(u_abc,theta,varargin)

    % Default settings
    Flag_ZeroSequence = false; % Default, no zero sequence
    Flag_PowerVariant = false; % Default, power invariant
    Flag_DirectMatrixCalculation = true; % Default, direct calculation
    
    % Check advanced settings
    for n = 1:length(varargin)
        if(strcmpi(varargin{n},'ZeroSequence'))
            Flag_ZeroSequence = varargin{n+1};
        elseif (strcmpi(varargin{n},'PowerVariant'))
            Flag_PowerInvariant = varargin{n+1};
        elseif (strcmpi(varargin{n},'DirectMatrixCalculation'))
            Flag_DirectMatrixCalculation = varargin{n+1};
        end
    end
    
 	% Check the form of input signal
    [r,c] = size(u_abc);
    if (r~=3) || (c~=1)
        error('Error: Input signal is not a 3*1 vector.')
    end
    
    % Calculate the transformation matrix directly
    if Flag_DirectMatrixCalculation
        % Get the transformation matrix
        if Flag_PowerVariant
         	T = 2/3*[cos(theta), cos(theta-2/3*pi), cos(theta-4/3*pi);
                     -sin(theta),-sin(theta-2/3*pi),-sin(theta-4/3*pi);
                     1/2,        1/2,               1/2];
        else
         	T = sqrt(2/3)*[cos(theta), cos(theta-2/3*pi), cos(theta-4/3*pi);
                           -sin(theta),-sin(theta-2/3*pi),-sin(theta-4/3*pi);
                           1/sqrt(2),  1/sqrt(2),         1/sqrt(2)];
        end
        
        % Check zero sequence
        if Flag_ZeroSequence
            T_abc2dq = T;
        else
            T_abc2dq = T(1:2,:);
        end
        
      	% Frame transformation
        u_dq = T_abc2dq * u_abc;
        
    % Call "abc2alphabeta" and "alphabeta2dq" to acheive two-step
    % transformation.
    else
        % Frame transformation
        u_alphabeta = SimplusGT.abc2alphabeta(u_abc,'PowerVariant',Flag_PowerVariant,'ZeroSequence',Flag_ZeroSequence);
        u_dq = SimplusGT.alphabeta2dq(u_alphabeta,theta);
    end

end