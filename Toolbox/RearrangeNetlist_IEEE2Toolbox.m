% This function transfers netlise of lines from IEEE form to toolbox form
% if required.

% Author(s): Yitong Li

function [UpdateLine] = RearrangeNetlist_IEEE2Toolbox(ListLine,ListLineIEEE)

if ListLineIEEE(1,1) == 1
    % "ListLineIEEE" is enabled, which over-writes "ListLine"
    ListLineIEEE = ListLineIEEE(4:end,:);
    
    % Organize mutual branch data
    FB_m = ListLineIEEE(:,1);
    TB_m = ListLineIEEE(:,2);
    Rbr_m = ListLineIEEE(:,3);
    Xbr_m = ListLineIEEE(:,4);
    Bbr_pi = ListLineIEEE(:,5);
    Gbr_pi = ListLineIEEE(:,6);
    Tbr_m = ListLineIEEE(:,7);
    
    % Check self branch error
    Dif_m = FB_m - TB_m;
    if ~isempty(find(~Dif_m, 1))
        error(['Error: Netlist IEEE can not have self branch.'])
    end
    
    % Number of bus
    N_Branch_m = length(FB_m);              % Number of mutual branch
    N_Bus_ = max(max(FB_m),max(TB_m));
    N_Branch_s = N_Bus_;                    % Number of self branch
    
    % Get the self branch data
    FB_s = transpose([1:N_Branch_s]);
    TB_s = FB_s;
    
    % Initialize R, X, B, G for self branch
    Rbr_s = zeros(N_Branch_s,1);
    Xbr_s = Rbr_s;
    Bbr_s = Rbr_s;
    Gbr_s = Rbr_s;
    Tbr_s = ones(N_Branch_s,1);
    
    for i = 1:N_Branch_m
        Bbr_s(FB_m(i)) = Bbr_s(FB_m(i))+Bbr_pi(i);
        Bbr_s(TB_m(i)) = Bbr_s(TB_m(i))+Bbr_pi(i);
     	Gbr_s(FB_m(i)) = Gbr_s(FB_m(i))+Gbr_pi(i);
        Gbr_s(TB_m(i)) = Gbr_s(TB_m(i))+Gbr_pi(i);
    end
    
    % Use IEEE form to over-write toolbox form 
    UpdateLine = ListLine;
else
    % Use the toolbox form
    UpdateLine = ListLine;
end

end