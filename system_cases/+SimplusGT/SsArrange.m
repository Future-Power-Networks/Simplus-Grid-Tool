% This function achieves the "arrange" operation of two state space
% systems.

% Author(s): Yitong Li

function Gss = SsArrange(GssCell)

[rmax,cmax] = size(GssCell);

for r = 1:rmax
    for c = 1:cmax
        if is_dss(GssCell{r,c})
            error(['Error: System 1 and/or 2 is not in ss form.']); 
        end
        GdssCell{r,c} = ss2dss(GssCell{r,c});
    end
end

Gdss = dss_Arrange(GdssCell);
Gss = dss2ss(Gdss);

end