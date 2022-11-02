load EigHz_1.mat
eig1=EigVecHz;
load EigHz_X1.mat
eigX1=EigVecHz;
load EigHz_2.mat
eig2=EigVecHz;
load EigHz_X2.mat
eigX2=EigVecHz;
load EigHz_3.mat
eig3=EigVecHz;
load EigHz_X3.mat
eigX3=EigVecHz;


figure(1)
clf
subplot(3,2,1)
scatter(real(eigX1),imag(eigX1),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eig1),imag(eig1),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Case-1 pole map');
axis([-40,10,-100,100]);
legend({'Without IBR2','With IBR2'})
subplot(3,2,2)
scatter(real(eigX1),imag(eigX1),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eig1),imag(eig1),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Case-1 pole map (zoomed)');
axis([-3,1,-100,100]);
legend({'Without IBR2','With IBR2'})

subplot(3,2,3)
scatter(real(eigX2),imag(eigX2),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eig2),imag(eig2),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Case-2 pole map');
axis([-40,10,-100,100]);
legend({'Without IBR2','With IBR2'})
subplot(3,2,4)
scatter(real(eigX2),imag(eigX2),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eig2),imag(eig2),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Case-2 pole map (zoomed)');
axis([-3,1,-100,100]);
legend({'Without IBR2','With IBR2'})

subplot(3,2,5)
scatter(real(eigX3),imag(eigX3),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eig3),imag(eig3),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Case-3 pole map');
axis([-40,10,-100,100]);
legend({'Without IBR2','With IBR2'})
subplot(3,2,6)
scatter(real(eigX3),imag(eigX3),'x','LineWidth',1.5); hold on; grid on;
scatter(real(eig3),imag(eig3),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Case-3 pole map (zoomed)');
axis([-3,1,-100,100]);
legend({'Without IBR2','With IBR2'})