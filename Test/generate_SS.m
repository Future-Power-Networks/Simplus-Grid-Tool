% This function is used to generate the state space matrix based on the
% input and output equations and state space vectors.

% State space model:
% x' = Ax + Bu
% y = Cx + Du

% u: input vector
% x: state vector
% y: output vector

% input_M: input equation from x, u to x', i.e. x' = input_M*[x;u];
% output_M: output equation from x, u to y, i.e. y = output_M*[x;u];

function [A, B, C, D] = generate_SS(input_M,output_M,u,x,y)

lu = length(u);
lx = length(x);
ly = length(y);

A = input_M(:,1:lx);
B = input_M(:,(lx+1):(lx+lu));

C = output_M(:,1:lx);
D = output_M(:,(lx+1):(lx+lu));

% for i = 1:lx;
%     for j = 1:lx
%        	A(i,j) = input_M(i,j);
%     end
%     for j = 1:lu
%         B(i,j) = input_M(i,j+lx);
%     end
% end
% 
% for m = 1:ly;
%     for n = 1:lx
%        	C(m,n) = output_M(m,n);
%     end
%     for n = 1:lu
%         D(m,n) = output_M(m,n+lx);
%     end
% end

end