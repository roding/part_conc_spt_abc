clear
clc
close all hidden

n1 = 10000;
n2 = 10000;

x1 = 3 + poissrnd(2, [n1, 1]);
x2 = x1%3 + poissrnd(2, [n2, 1]);

y1 = lognrnd(2, 0.01, [n1, 1]);
y2 = y1%lognrnd(2, 0.01, [n2, 1]);

d = 0;

tic
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
toc

% tic
% X = [x1 -ones(size(x1)); x2 +ones(size(x2))];
% [~, ind] = sort(X(:, 1), 'ascend');
% X(:, 1) = X(ind, 1);
% X(:, 2) = X(ind, 2);
% p = cumsum(X(:,2)); plot(p)
% sum(p.^2)
% mean(x1)
% mean(x2)
tic
x1 = sort(x1);
x2 = sort(x2);
if n1 <= n2
    p = sum((x1 - x2(1:n1)).^2) + sum(x2(n1+1:end).^2);
else
    p = sum((x1(1:n2) - x2).^2) + sum(x1(n2+1:end).^2);
end
p
y1 = sort(y1);
y2 = sort(y2);
if n1 <= n2
    p = sum((y1 - y2(1:n1)).^2) + sum(y2(n1+1:end).^2);
else
    p = sum((y1(1:n2) - y2).^2) + sum(y1(n2+1:end).^2);
end
p

toc
