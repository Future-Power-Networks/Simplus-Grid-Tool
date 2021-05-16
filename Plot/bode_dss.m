% This function plots bode diagram for a descriptor state space system.

% Authors(s): Yitong Li

%%
% Notes:
% 
% The system can be either in state space form or descriptor state space
% form; The system can be either in real state space form or complex state
% space form.
%
% Input of this function:
% X        - target system in descriptor state space form
% jw       - frequency range
% varargin - Variable length input argument list, for advanced settings

%%
function Xw = bode_dss(X,jw,varargin)

    [PhaseOn,~] = LoadVar(1,'PhaseOn',varargin);    % Default, 1
    [PlotOn,~]  = LoadVar(1,'PlotOn',varargin);   	% Default, 1

    [M,N] = size(X);

    % Calculate the system gain at the given frequency range        
    % Initialize Xw
    Xw = zeros(M,N,length(jw));
    % Calculate Xw
    for n = 1:length(jw)
        Xw(:,:,n) = dss2value(X,jw(n));
    end
    % Check inf and NaN
    Xw_nan = isnan(Xw);
    Xw_inf = isinf(Xw);
    if ~isempty(find(Xw_nan,1)); error(['Error: system gain is NaN']); end
	if ~isempty(find(Xw_inf,1)); error(['Error: system gain is inf']); end
    
    % Plot the bode diagram
    if PlotOn == 1
        plotc(Xw,imag(jw)/(2*pi),'PhaseOn',PhaseOn,varargin);
    end
    
end
    
    


