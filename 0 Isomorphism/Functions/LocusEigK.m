scale = logspace(-2,0,10)
% scale = linspace(1,0.01,10);

Hinv_ = Hinv;
for i = 1:length(scale)
    Hinv_(1) = Hinv(1)*scale(i);
    KH_ = Hinv_*K
    xi_{i} = eig(KH_);
end

figure_n = 9002;
figure(figure_n)
for i = 1:length(scale)
scatter(real(xi_{i}),imag(xi_{i}),'x','LineWidth',1.5); hold on; grid on;
end