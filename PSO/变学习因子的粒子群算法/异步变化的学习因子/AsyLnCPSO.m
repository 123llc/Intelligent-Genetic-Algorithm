%% 变学习因子的粒子群算法-异步变化的学习因子
% 粒子数目：N
% 学习因子1的最大值：c1max
% 学习因子1的最小值：c1min
% 学习因子2的最大值：c2max
% 学习因子2的最小值：c2min
% 惯性权重 w
% 最大迭代次数 M 
% 自变量的个数 D
% 目标函数取值最小时的自变量值 xm
% 目标函数的最小值 fv
% 调用格式：[xm,fv] = AsyLnCPSO(fitness,N,c1max,clmin,c2max,c2min,w,M,D) 

function [xm,fv] = AsyLnCPSO(fitness,N,c1max,c1min,c2max,c2min,w,M,D)

format long; % 设置输出格式
for i=1:N
    for j=1:D
        x(i,j) = randn; %初始化速度和位置
        v(i,j) = randn;
    end
end

for i=1:N
    p(i) = fitness(x(i,:)); % p(i) 是产生的适应度的值
    y(i,:) = x(i,:); % y(i,:)，x(i,:) 保存的位置的值
end

pg = x(N,:);  % pg为全局最优
for i=1:(N-1)
    if fitness(x(i,:))<fitness(pg)
        pg = x(i,:);
    end
end

for t = 1:M
    c1 = c1max - (c1max - c1min)*t/M; %线性变化的学习因子1
    c2 = c2max - (c2max - c2min)*t/M; %线性变化的学校因子2
    for i = 1:N
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        if fitness(x(i,:)) < p(i)
            p(i) = fitness(x(i,:));
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
        