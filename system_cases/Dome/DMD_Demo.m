clear all;
close all;
clc;

%% generate grid geometry
xi = linspace(-10, 10, 400);
t = linspace(0, 4*pi, 200);
dt = t(2) - t(1);

[xgrid, tgrid] = meshgrid(xi, t);

%% create two spatio-temporal patterns
f1 = sech(xgrid + 3) .* (1*exp(1i*2.3*tgrid));
f2 = (sech(xgrid).*tanh(xgrid)).*(2*exp(1i*2.8*tgrid));

f = f1 + f2;

[u, s, v] = svd(f.');

%% plot
figure(1);
subplot(2,2,1); surfl(xgrid, tgrid, real(f1));  shading interp; colormap gray;  title('f1(x, t)');
subplot(2,2,2); surfl(xgrid, tgrid, real(f2));  shading interp; colormap gray;  title('f2(x, t)');
subplot(2,2,3); surfl(xgrid, tgrid, real(f));  shading interp; colormap gray; title('f(x, t) = f1(x, t) + f2(x, t)');

figure(2);
plot(diag(s) / sum(diag(s)), 'ro'); title('SVD: low rank property (rank = 2, two modes)');

figure(3);
subplot(2,1,1); plot(real(u(:, 1:3))); 
legend('1st mode of basis u (left singular vectors)', ...
        '2nd mode of basis u (left singular vectors)', ...
        '3rd mode; numerical round off (junk basis)  '); 
title('Modes for basis of column space of '' f ''');

subplot(2,1,2); plot(real(v(:, 1:3)));
legend('1st mode of basis v (right singular vectors)', ...
        '2nd mode of basis v (right singular vectors)', ...
        '3rd mode; numerical round off (junk basis)  '); 
title('Modes for basis of row space of '' f ''');

%% dynamic mode decomposition (DMD)
X = f.';    % in C^(spatio, temporal)

X1 = X(:, 1:end-1);
X2 = X(:, 2:end);

%% STEP 1: singular value decomposition (SVD)
r = 2;      % rank-r truncation
[U, S, V] = svd(X1, 'econ');

Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);

%% STEP 2: low-rank subspace matrix 
%         (similarity transform, least-square fit matrix, low-rank subspace matrix)
Atilde = Ur'*X2*Vr*Sr^(-1);

%% STEP 3: eigen decomposition
% W: eigen vectors
% D: eigen values
[W, D] = eig(Atilde);

%% STEP 4: real space DMD mode
Phi = X2*Vr*Sr^(-1)*W;   % DMD modes

lamdba = diag(D);       % eigen value
omega = log(lamdba)/dt; % log of eigen value

figure(4); 
subplot(2,1,1); plot(real(u(:, 1:2)));
legend('1st mode of SVD', ...
        '2nd mode of SVD'); 
subplot(2,1,2); plot(real(Phi));
legend('1st mode of DMD', ...
        '2nd mode of DMD'); 

%% STEP 5: reconstruct the signal
x1 = X(:, 1);       % time = 0
b = pinv(Phi)*x1;   % initial value; \: pseudo inverse

t_dyn = zeros(r, length(t));

for i = 1:length(t)
   t_dyn(:, i) = (b.*exp(omega*t(i))); 
end

f_dmd = Phi*t_dyn;

figure(1);
subplot(2,2,1); surfl(xgrid, tgrid, real(f1));  shading interp; colormap gray;  title('f1(x, t)');
subplot(2,2,2); surfl(xgrid, tgrid, real(f2));  shading interp; colormap gray;  title('f2(x, t)');
subplot(2,2,3); surfl(xgrid, tgrid, real(f));  shading interp; colormap gray; title('f(x, t) = f1(x, t) + f2(x, t)');
subplot(2,2,4); surfl(xgrid, tgrid, real(f_dmd).');  shading interp; colormap gray; title('reconstruction of f(x, t) by DMD');

%% additional informations
% predict furture state using time dynamics
t_ext = linspace(0, 8*pi, 400);

[xgrid_ext, tgrid_ext] = meshgrid(xi, t_ext);

t_ext_dyn = zeros(r, length(t_ext));

for i = 1:length(t_ext)
   t_ext_dyn(:, i) = (b.*exp(omega*t_ext(i))); 
end

f_dmd_ext = Phi*t_ext_dyn;

figure(4);
subplot(2,2,1); surfl(xgrid, tgrid, real(f));  shading interp; colormap gray; 
xlabel('spatial-axis'); ylabel('temporal-axis'); title('f(x, t) during t = [0, 4*pi]');
subplot(2,2,2); surfl(xgrid, tgrid, real(f_dmd).');  shading interp; colormap gray; 
xlabel('spatial-axis'); ylabel('temporal-axis'); title('DMD reconstruction of f(x, t) during t = [0, 4*pi]');
subplot(2,2,[3,4]); surfl(xgrid_ext, tgrid_ext, real(f_dmd_ext).');  shading interp; colormap gray; 
xlabel('spatial-axis'); ylabel('temporal-axis'); title('DMD prediction of f(x, t) during t = [0, 8*pi]');

% If eigen value: lambda or omega has tiny real part > 0,
% the output of system function which is spaned 
% by eigen vectors with eigen values goes to the infinity.
% It is one of the limitation of DMD method.

format long e;
disp(omega);    % eigen values exist on the imaginary part.
disp('real part: numerical round off');
disp('imag part: eigen values');