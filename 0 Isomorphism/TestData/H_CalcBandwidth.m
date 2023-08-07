% This function is used to calculate the bandwidth with a given lambda

function wb = H_CalcBandwidth(lambda_w)

f0 = 60;
w0 = f0*2*pi;

wb_vec = linspace(0.2,3,10000)*w0;

for k = 1:length(wb_vec)
F = (1i*w0 - lambda_w)/(1i*(w0 + wb_vec(k)) - lambda_w);   % Low pass filter
Mag_F = abs(F);
Angle_F = angle(F);
    if (Mag_F <= 1/sqrt(2)) || (Angle_F >= pi/4) || (Angle_F<=-pi/2)
        wb = wb_vec(k);
        fb = wb/2/pi;
        break;
    end
end

end