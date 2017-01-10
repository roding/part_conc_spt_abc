clear
clc
close all hidden

number_of_columns = 5;
precision_in_bytes = 8;
    
folder = 'res_20170109';
files = dir([folder '/' '*.dat']);
number_of_files = numel(files);

number_of_files = 100;

mu = [];
sigma = [];
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

%     mu{current_file} = data(:, 1);
%     sigma{current_file} = data(:, 2);
%     c{current_file} = data(:, 3);
%     az{current_file} = data(:, 4);
%     dist{current_file}  = data(:, 5);
    
%     mu = [mu ; data(:, 1)];
%     sigma =  [sigma ; data(:, 2)];
%     c = [c ; data(:, 3)];
%     az = [az ; data(:, 4)];
%     dist = [dist ; data(:, 5)];
    
    mu = cat(1, mu, data(:, 1));
    sigma = cat(1, sigma, data(:, 2)); 
    c = cat(1, c, data(:, 3));
    az = cat(1, az, data(:, 4));
    dist = cat(1, dist, data(:, 5));
end

%% Try inference.
p = 0.00005;
epsilon = prctile(dist, p * 100);
index = dist <= epsilon;

mu = mu(index);
sigma = sigma(index);
c = c(index);
az = az(index);

figure, hold on, hist(mu, 100), title('mu'), hold off
figure, hold on, hist(sigma, 100), title('sigma'), hold off
figure, hold on, hist(c, 100), title('c'), hold off
figure, hold on, hist(az, 100), title('az'), hold off