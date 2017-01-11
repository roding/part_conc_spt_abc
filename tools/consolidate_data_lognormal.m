clear
clc
close all hidden

number_of_columns = 5;
precision_in_bytes = 8;
    
folder = '../test/res';
files = dir([folder '/' '*.dat']);
number_of_files = numel(files);


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
p = 0.005;
epsilon = prctile(dist, p * 100);
index = dist <= epsilon;

mu = mu(index);
sigma = sigma(index);
c = c(index);
az = az(index);

m = exp(mu + 0.5*sigma.^2);
s = m .* sqrt(exp(sigma.^2) - 1);

figure, hold on, hist(mu, 50), title('mu'), hold off
figure, hold on, hist(sigma, 50), title('sigma'), hold off
figure, hold on, hist(m, 50), title('m'), hold off
figure, hold on, hist(s, 50), title('s'), hold off
figure, hold on, hist(c, 50, 100), title('c'), hold off
figure, hold on, hist(az, 50, 100), title('az'), hold off