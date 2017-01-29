clear
clc
close all hidden

%% Read data.
file_name = '../test/abc_pmc_sample_lognormal_1.dat';
file_info = dir(file_name);
file_size = file_info.bytes;
number_of_columns = 4;
precision_in_bytes = 8;
number_of_abc_samples = file_size / number_of_columns / precision_in_bytes;

file_id = fopen(file_name);
data = fread(file_id, [number_of_columns, number_of_abc_samples], 'float64');
fclose(file_id);

data = data';

m = data(:, 1);
s = data(:, 2);
c = data(:, 3);
az = data(:, 4);
% dist  = data(:, 5);

clear data

%% Try inference.
% p = 0.005;
% epsilon = prctile(dist, p * 100);
% index = dist <= epsilon;
% 
% mu = mu(index);
% sigma = sigma(index);
% c = c(index);
% az = az(index);
% 
% m = exp(mu + 0.5*sigma.^2);
% s = m .* sqrt(exp(sigma.^2) - 1);

% figure, hold on, hist(mu), title('mu'), hold off
% figure, hold on, hist(sigma), title('sigma'), hold off
figure, hold on, hist(m), title('m'), hold off
figure, hold on, hist(s), title('s'), hold off
figure, hold on, hist(c), title('c'), hold off
figure, hold on, hist(az), title('az'), hold off