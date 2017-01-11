clear
clc
close all hidden

%% Read data.
file_name = 'simulated_system.dat';

data = dlmread(file_name, ',');

% data = data';

K = data(:, 1);
DE = data(:, 2);

clear data

%% Analyze.

hist(DE, 1000)

mean(DE)