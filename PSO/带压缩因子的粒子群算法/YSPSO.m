%% 带压缩因子的粒子群算法
% 粒子数目：N
% 学习因子 c1,c2
% 惯性权重 w
% 最大迭代次数 M 
% 自变量的个数 D
% 目标函数取值最小时的自变量值 xm
% 目标函数的最小值 fv

% 调用格式：[xm,fv] = YSPSO(fitness,N,c1,c2,w,M,D) 

function [xm,fv] = YSPSO(fitness,N,c1,c2,M,D)

phi = c1+c2;
if phi <= 4
    disp('c1与c2的和必须大于4！');
    xm = NaN;
    fv = NaN;
    return;
end
format long; % 设置输出格式
for i=1:N
    for j=1:D
        x(i,j) = randn;
        v(i,j) = randn;
    end
end

for i=1:N
    p(i) = fitness(x(i,:));
    y(i,:) = x(i,:);
end
pg = x(N,:);
for i=1:(N-1)
    if fitness(x(i,:))<fitness(pg)
        pg = x(i,:);
    end
end
for t = 1:M
    for i = 1:N
        ksi = 2/abs(2-phi-sqrt(phi^2 - 4*phi)); % 压缩因子
        v(i,:) = v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        v(i,:) = ksi*v(i,:);
        x(i,:) = x(i,:)+v(i,:);
        if fitness(x(i,:)) < p(i)
            p(i) =fitness(x(i,:));
            y(i,:) = x(i,:);
        end
        % 每次循环都要比较
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
end
% '代表转置的意思
xm = pg';
fv = fitness(pg);
        