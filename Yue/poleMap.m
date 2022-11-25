load EigHz_C1_200.mat
eigC1_200=EigVecHz;
load EigHz_C2_200.mat
eigC2_200=EigVecHz;


figure(1)
clf
scatter(real(eigC1_200),imag(eigC1_200),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eigC2_200),imag(eigC2_200),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('pole map');
axis([-40,10,-100,100]);
legend({'case1','case2'})


load EigHz_C1_200.mat
eigC1_200=EigVecHz;
load EigHz_C3_200.mat
eigC3_200=EigVecHz;


figure(2)
clf
scatter(real(eigC1_200),imag(eigC1_200),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eigC3_200),imag(eigC3_200),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('pole map');
axis([-40,10,-100,100]);
legend({'case1','case3'})