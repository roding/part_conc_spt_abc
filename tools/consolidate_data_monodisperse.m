clear
clc
close all hidden

number_of_columns = 4;
precision_in_bytes = 8;
    
folder = 'res_20170111';
files = dir([folder '/' '*.dat']);
number_of_files = numel(files);

D = [];
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
    
    D = cat(1, D, data(:, 1));
    c = cat(1, c, data(:, 2));
    az = cat(1, az, data(:, 3));
    dist = cat(1, dist, data(:, 4));
end

%% Try inference.
p = 0.001;
epsilon = prctile(dist, p * 100);
index = dist <= epsilon;

D = D(index);
c = c(index);
az = az(index);

figure, hold on, hist(D, 100), title('D'), hold off
figure, hold on, hist(c, 100), title('c'), hold off
figure, hold on, hist(az, 100), title('az'), hold off