function thetaNew = DisableMod(theta)
    len = length(theta);
    theta(len+1) = theta(len);
    for i = 1:len
        if theta(i+1) - theta(i) > pi
            theta(i+1:end) = theta(i+1:end) - 2*pi;
        elseif theta(i+1) - theta(i) < -pi
            theta(i+1:end) = theta(i+1:end) + 2*pi;
        end
    end
    thetaNew = theta(1:len);
end