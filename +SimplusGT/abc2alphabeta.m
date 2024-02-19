% This function transforms a time-domain signal from natural frame (abc
% frame) to statationary frame (alpha/beta/0 frame).

% Author(s): Yitong Li

function u_alphabeta0 = abc2alphabeta(u_abc,varargin)

    % Default settings
    Flag_ZeroSequence = false;  % No zero sequence
    Flag_PowerVariant = false;  % Power invariant

    % Check advanced settings
    for n = 1:length(varargin)
        if(strcmpi(varargin{n},'ZeroSequence'))
            Flag_ZeroSequence = varargin{n+1};
        elseif (strcmpi(varargin{n},'PowerVariant'))
            Flag_PowerVariant = varargin{n+1};
        end
    end
    
    % Check the form of input signal
    [r,c] = size(u_abc);
    if (r~=3) || (c~=1)
        error('Error: Input signal is not a 3*1 vector.')
    end

    % Get the Clarke transformation matrix
	if Flag_PowerVariant
        % Power variant
        T_Clarke = 2/3*[1,   -1/2,      -1/2;
                        0,   sqrt(3)/2, -sqrt(3)/2;
                        1/2, 1/2,       1/2];
    else
        % Power invariant
    	T_Clarke = sqrt(2/3)*[1,         -1/2,      -1/2;
                              0,         sqrt(3)/2, -sqrt(3)/2;
                              1/sqrt(2), 1/sqrt(2), 1/sqrt(2)];
    end

    if Flag_ZeroSequence
        T_abc2alphabeta = T_Clarke;
    else
        T_abc2alphabeta = T_Clarke(1:2,:);
    end
                  
    % Frame transformation
    u_alphabeta0 = T_abc2alphabeta*u_abc;

end