clear
% clc
close all hidden

%% Read data.
file_name = '../test/simulated_system.dat';

data = dlmread(file_name, ',');

K = data(:, 1);
DE = data(:, 2);

clear data

%% Analyze.

figure
hist(DE, 1000);

figure
hist(K, 1:100);

mean(DE)
std(DE)