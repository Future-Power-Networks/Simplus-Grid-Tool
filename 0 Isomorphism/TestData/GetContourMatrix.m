% This function gets the contour matrix of the channel low pass filter of a
% power grid, so that it can be plotted directly by using the matlab
% "contour" function.

% Author(s): Yitong Li

function ContourMatrix = GetContourMatrix(fb,x_vec,y_vec)

%% Funfamentals
f0 = 60;                % We consider the US network here, so f0 should be 60Hz.
w0 = 2*pi*f0;

wb = 2*pi*fb;

x = sym('x','real');    % Real part of lambda
y = sym('y','real');    % Imaginary part of lambda

lambda_f = x + 1i*y;
lambda_w = 2*pi*lambda_f;

F = (1i*w0 - lambda_w)/(1i*(w0 + wb) - lambda_w);   % Low pass filter

%% Get the function
Real_F = real(F);
Imag_F = imag(F);

Real_Fun = matlabFunction(Real_F);
Imag_Fun = matlabFunction(Imag_F);

Mag_F_sqr = Real_F^2 + Imag_F^2;
Mag_F_sqr = Mag_F_sqr - 1/2;        % The magnitude is compared to sqrt(2)
Mag_Fun = matlabFunction(Mag_F_sqr);

Ang_F = atan2(Imag_F,Real_F);
Ang_Fun = matlabFunction(Ang_F);
Ang_F1 = Ang_F - pi/2/2;            % The phase angle is compared to +-45 degree
Ang_F2 = -pi/2/2 - Ang_F;
Ang_Fun1 = matlabFunction(Ang_F1);
Ang_Fun2 = matlabFunction(Ang_F2);

%% Get the matrix
xlength = length(x_vec);
ylength = length(y_vec);

% Get F_Sign matrix
% F_Sign is the matrix that classifies the access region (labelled by 1)
% and forbidden region (labelled by -1).
for n = 1:ylength
    for m = 1:xlength
        % Get the x,y with a very small perturbation, to avoid the NaN for
        % the case such as 0/0.
        x_ = x_vec(m) + eps(x_vec(m));
        y_ = y_vec(n) + eps(y_vec(n));
        
        % Calculate magnitude
        Mag_F_Value(n,m) = Mag_Fun(x_,y_);
        Mag_F_Sign(n,m) = sign(Mag_F_Value(n,m));
        
        % Calculate phase angle
        Ang_F_Value1(n,m) = Ang_Fun1(x_,y_);
        Ang_F_Value2(n,m) = Ang_Fun2(x_,y_);
        if (Ang_F_Value1(n,m)>0) || (Ang_F_Value2(n,m)>0)
            Ang_F_Sign(n,m) = -1;
        else
            Ang_F_Sign(n,m) = 1;
        end
        
        % Calculate the whole forbidden area
        % Notes: -1 is the forbidden region, i.e., the area of lambda that
        % does not satisfy the requirement of bandwidth wb.
        if (Mag_F_Sign(n,m)<0) || (Ang_F_Sign(n,m)<0)
            F_Sign(n,m) = -1;
        else
            F_Sign(n,m) = 1;
        end
    end
end

% Get F_Plot
% F_Plot is the matrix that can be plotted by matlab contour function. The
% logic of finding F_Plot is to find the boundary that connects 1 and -1 in
% F_Sign matrix.
F_Plot = zeros(ylength,xlength);
for n = 1:ylength
    for m = 1:xlength
        if F_Sign(n,m) == 1
            Position = [n,m];
            value = 1;
            mode = 2;
            Check = CheckNeighbor(F_Sign,Position,value,mode);
            if Check ~= 1
                F_Plot(n,m) = 1;
            end
        end
    end
end

ContourMatrix = F_Plot;

end