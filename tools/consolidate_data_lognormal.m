clear
clc
close all hidden

number_of_columns = 5;
precision_in_bytes = 8;
    
folder = '../test/res_large_system_m_s_20170113';
files = dir([folder '/' '*.dat']);
number_of_files = numel(files);


m = [];
s = [];
c = [];
az = [];
dist = [];
    
for current_file = 1:number_of_files
    disp(current_file)
    file_path =  [folder '/' files(current_file).name];
    file_info = dir(file_path);
    file_size = file_info.bytes;

    number_of_abc_samples = file_size / number_of_columns / precision_in_bytes;

    file_id = fopen(file_path);
    data = fread(file_id, [number_of_columns, number_of_abc_samples], 'float64');
    fclose(file_id);

    data = data';

    m = cat(1, m, data(:, 1));
    s = cat(1, s, data(:, 2)); 
    c = cat(1, c, data(:, 3));
    az = cat(1, az, data(:, 4));
    dist = cat(1, dist, data(:, 5));
end

%% Try inference.
p = 0.001;
epsilon = prctile(dist, p * 100);
index = dist <= epsilon;

m = m(index);
s = s(index);
c = c(index);
az = az(index);

figure, hold on, hist(m, 50), title('m'), hold off
figure, hold on, hist(s, 50), title('s'), hold off
figure, hold on, hist(c, 50, 100), title('c'), hold off
figure, hold on, hist(az, 50, 100), title('az'), hold off