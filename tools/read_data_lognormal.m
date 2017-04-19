clear
clc
close all hidden

file_name = '../src/output_lognormal/res_lognormal_323.dat';
file_info = dir(file_name);
file_size = file_info.bytes;
number_of_columns = 6;
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
dist  = data(:, 5);
w  = data(:, 6);

clear data


figure, hold on, hist(m), title('m'), hold off
figure, hold on, hist(s), title('s'), hold off
figure, hold on, hist(c), title('c'), hold off
figure, hold on, hist(az), title('az'), hold off
%figure, hold on, scatter(m, w)
figure, hold on, scatter(s, log10(dist))
mean(m)
mean(s)
mean(c)
mean(az)
max(log10(dist))