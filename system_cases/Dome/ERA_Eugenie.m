V2 = timeseries2timetable(out.V_log{1}.Values);
V2 = V2(5988000:6036000,:);%6120000
V2.Properties.VariableNames(1) = "data";
V2.data = 1000000*(V2.data - V2.data(1));
V2 = V2(12000:end,:);

figure
plot(V2.Time,V2.data)
hold on

f_sample_new = 1600;
V2_rt = retime(V2,'regular','SampleRate',f_sample_new);
V2_rt.data(1) = 0;
plot(V2_rt.Time,V2_rt.data)
%%
sys = era(V2_rt,30);
e = eig(sys.A);
%%
N = length(V2_rt.data)-1;
L = 180;
H1 = hankel(V2_rt.data(1:N-L+1),V2_rt.data(N-L+1:N));
H2 = hankel(V2_rt.data(2:N-L+2),V2_rt.data(N-L+2:N+1));
[U1,S1,V1] = svd(H1);
[U2,S2,V2] = svd(H2);
new_order = 5;
U1_reduced = U1(:,1:new_order);
S1_reduced = S1(1:new_order,1:new_order);
V1_reduced = V1(:,1:new_order);
U2_reduced = U2(:,1:new_order);
S2_reduced = S2(1:new_order,1:new_order);
V2_reduced = V2(:,1:new_order);
H2_reduced = U2_reduced*S2_reduced*V2_reduced.';
% obs_mtx = U1_reduced*S1_reduced^(1/2);
% ctrl_mtx = S1_reduced^(1/2)*V1_reduced.';
A = S1_reduced^(-1/2)*U1_reduced.'*H2_reduced*V1_reduced*S1_reduced^(-1/2);
e = eig(A);
%%
figure;
subplot(1,2,1)
scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
scatter(real(e),imag(e),'o','LineWidth',1.5);
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Global pole map');

subplot(1,2,2)
scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
scatter(real(e),imag(e),'o','LineWidth',1.5);
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Zoomed pole map');
axis([-80,20,-150,150]);
plot([-80,0], [-80,0]*10, '--k','LineWidth',2,'Color','blue')
plot([-80,0], [80,0]*10, '--k','LineWidth',2,'Color','blue')
legend('mode','ERA','10% damping line')