clear
clc
close all hidden

file_name = '../src/output_monodisperse/res_monodisperse_424.dat';
file_info = dir(file_name);
file_size = file_info.bytes;
number_of_columns = 5;
precision_in_bytes = 8;
number_of_abc_samples = file_size / number_of_columns / precision_in_bytes;

file_id = fopen(file_name);
data = fread(file_id, [number_of_columns, number_of_abc_samples], 'float64');
fclose(file_id);

data = data';

D = data(:, 1);
c = data(:, 2);
az = data(:, 3);
dist = data(:, 4);
w = data(:, 5);

clear data


figure, hold on, hist(D), title('D'), hold off
figure, hold on, hist(c), title('c'), hold off
figure, hold on, hist(az), title('az'), hold off
mean(D)
mean(c)
mean(az)
max(log10(dist))