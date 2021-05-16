% This function converts the a complex value from polar form to rectangular form

% Author(s): Yunjie Gu， Yitong Li

function rect = pol2rect(r,o)   % r = magnitude, o = angle in radians.
rect = r.*cos(o) + j*r.*sin(o);   % rect = real + j*imag
end