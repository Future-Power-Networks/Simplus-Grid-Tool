% Tihs function transforms a time-domain signal from stationary frame
% (alpha/beta/0 frame) to natural frame (abc frame).

% Author(s): Yitong Li

function u_abc = alphabeta2abc(u_alphabeta,varargin)

    % Default settings
    Flag_PowerVariant = false; % Default, power invariant
    Flag_DirectMatrixCalculation = true; % Default, true
    
    % Check advanced settings
    for n = 1:length(varargin)
        if (strcmpi(varargin{n},'PowerVariant'))
            Flag_PowerVariant = varargin{n+1};
        elseif (strcmpi(varargin{n},'DirectMatrixCalculation'))
            Flag_DirectMatrixCalculation = varargin{n+1};
        end
    end
    
    % Check the form of input signal
    [r,c] = size(u_alphabeta);
    if (r==2) && (c==1)
        Flag_ZeroSequence = 0;
    elseif (r==3) && (c==1)
        Flag_ZeroSequence = 1;
    else
        error('Error: Input signal is not a 2*1 or 3*1 vector.')
    end

  	% Get the transformation matrix value directly. 
  	% This methed could eliminate the mathmatical calculation errors of
 	% Matlab, and therefore is recommended.
    if Flag_DirectMatrixCalculation

        if Flag_PowerVariant
            T_inv = [1,   0,           1;
                     -1/2, sqrt(3)/2,  1;
                     -1/2, -sqrt(3)/2, 1];
        else
            T_inv = sqrt(2/3)*[1,    0,          1/sqrt(2);
                               -1/2, sqrt(3)/2,  1/sqrt(2);
                               -1/2, -sqrt(3)/2, 1/sqrt(2)];
        end
        
 	% Calculate the transformation matrix by inverting the Clarke
 	% transformation.
    else
            % Get the Clarke transformation matrix
            if Flag_PowerVariant
                K = 2/3;        % Power variant
            else
                K = sqrt(2/3);  % Power invariant
            end
            T = [1,         -1/2,      -1/2;
                 0,         sqrt(3)/2, -sqrt(3)/2;
                 1/sqrt(2), 1/sqrt(2), 1/sqrt(2)];
            T_Clarke = K*T;

            % Get the inverse Clarke transformation matrix
            T_inv = inv(T_Clarke);
    end
    
    % Check zero sequence
  	if Flag_ZeroSequence
       	T_alphabeta2abc = T_inv;
    else
      	T_alphabeta2abc = T_inv(:,1:2);  % No zero sequence
  	end
                  
    % Frame transformation
    u_abc = T_alphabeta2abc*u_alphabeta;

end