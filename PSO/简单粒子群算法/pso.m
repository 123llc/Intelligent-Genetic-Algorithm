%% 基本粒子群算法
% 粒子数目：N
% 学习因子 c1,c2
% 惯性权重 w
% 最大迭代次数 M 
% 自变量的个数 D
% 目标函数取值最小时的自变量值 xm
% 目标函数的最小值 fv

% 调用格式：[xm,fv] = pso(fitness,N,c1,c2,w,M,D) 

function [xm,fv] = pso(fitness,N,c1,c2,w,M,D)

format long;
% 初始化速度和位置
for i = 1:N
    for j = 1:D
        x(i,j) = randn;
        v(i,j) = randn;
    end
end
% p(i) 是产生的适应度的值
% y(i,:)，x(i,:) 保存的位置的值
for i = 1:N
    p(i) = fitness(x(i,:)); % 适应度值
    y(i,:) = x(i,:);
end

pg = x(N,:);  % pg为全局最优
for i = 1:(N-1)
    if fitness(x(i,:) < fitness(pg))
        pg = x(i,:);
    end
end

% 速度、位置更新
for t = 1:M
    for i = 1:N
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        % 适应度值的比较
        if fitness(x(i,:)) < p(i)
            y(i,:)=x(i,:);
        end
        % 位置信息的比较
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
    % 最后的位置
    Pbest(t) = fitness(pg);
end
xm = pg';
fv = fitness(pg);
        
