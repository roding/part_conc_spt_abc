clear
clc
close all hidden

n1 = 100;
n2 = 100;

x1 = 3 + poissrnd(2, [n1, 1]);
x2 = 3 + poissrnd(2, [n2, 1]);

y1 = lognrnd(2, 0.01, [n1, 1]);
y2 = lognrnd(2, 0.01, [n2, 1]);

d = 0;

for i = 1:n1
    e1 = sum( x1 <= x1(i) & y1 <= y1(i) );
    e2 = sum( x2 <= x1(i) & y2 <= y1(i) );
    
    d = d + (e1-e2)^2;
end

for i = 1:n2
    e1 = sum( x1 <= x2(i) & y1 <= y2(i) );
    e2 = sum( x2 <= x2(i) & y2 <= y2(i) );
    
    d = d + (e1-e2)^2;
end

d