num = [1 2]; % Numerator coefficients
den = [1 4 4]; % Denominator coefficients
G = tf(num, den); % Create the transfer function

% Reduce the order using Pade approximation
order_reduction = 2;
G_pade = pade(G, order_reduction);
