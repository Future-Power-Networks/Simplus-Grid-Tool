rng(0);
% make dummy data
ncases = 30;
N = 1000;
data = NaN(ncases, N);
for k = 1:ncases
    % generate a dummy random-normal distribution
    % with random mean and random standard deviation
    data(k, :) = randn(N, 1);
end
% generate data for histograms: 1 histogram per column
edges = [-5:0.5:5]; % bin edges
counts = histc(data, edges, 2); % specify dim 2 to act column-wise
% plot results
hf = figure;
ha = axes;
hb = bar3(edges, counts.'); % note the transpose to get the colors right
xlabel('case number')
ylabel('bins');
zlabel('count');