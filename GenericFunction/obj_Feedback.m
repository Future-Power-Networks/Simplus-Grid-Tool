% This function achieves "dss_Feedback" when inputs are two matlab objects.
% Strings of state, input, and output are also re-arranged.

% Author(s): Yitong Li

%%
function Gobj = obj_Feedback(Gobj1,Gobj2,feedin,feedout,sign)

% Load object date
[~,G1] = Gobj1.ReadDSS(Gobj1);
[~,G2] = Gobj2.ReadDSS(Gobj2);
[StateStr1,InputStr1,OutputStr1] = Gobj1.ReadString(Gobj1);
[StateStr2,InputStr2,OutputStr2] = Gobj1.ReadString(Gobj2);

% Default
if nargin == 2
    [ly1,lu1] = size(G1.D);
    feedin = [1:lu1]; 
    feedout = [1:ly1]
end
if nargin == 3
    error(['Error: feedout is missing'])
end
if nargin <= 5
    sign = -1;      % Default, negative feedback
end

% Feedback the model
G = dss_Feedback(G1,G2,feedin,feedout,sign);

% Connect string
% The relationship is determined inside "dss_Feedback"
StateStr = [StateStr1,StateStr2];
InputStr = InputStr1;
OutputStr = OutputStr1;

% Create a new object
Gobj = Class_Model_DSS;
Gobj.LoadModel(Gobj,G);
Gobj.WriteString(Gobj,StateStr,InputStr,OutputStr);
obj_CheckDim(Gobj);

end