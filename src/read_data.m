clear
clc
close all hidden

%% Read data.
file_name = 'abc_sample.dat';
file_info = dir(file_name);
file_size = file_info.bytes;
number_of_columns = 5;
precision_in_bytes = 8;
number_of_abc_samples = file_size / number_of_columns / precision_in_bytes;

file_id = fopen(file_name);
data = fread(file_id, [number_of_columns, number_of_abc_samples], 'float64');
fclose(file_id);

data = data';

mu = data(:, 1);
sigma = data(:, 2);
c = data(:, 3);
az = data(:, 4);
dist  = data(:, 5);

clear data

%% Try inference.
p = 0.0001;
epsilon = prctile(dist, p * 100);
index = dist <= epsilon;

mu = mu(index);
sigma = sigma(index);
c = c(index);
az = az(index);

figure, hold on, hist(mu), title('mu'), hold off
figure, hold on, hist(sigma), title('sigma'), hold off
figure, hold on, hist(c), title('c'), hold off
figure, hold on, hist(az), title('az'), hold off