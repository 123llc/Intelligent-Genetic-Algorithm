%% 变学习因子的粒子群算法-同步变化的学习因子
% 粒子数目：N
% 学习因子的最大值：cmax
% 学习因子的最小值：cmin
% 随机权重的方差
% 最大迭代次数 M 
% 自变量的个数 D
% 目标函数取值最小时的自变量值 xm
% 目标函数的最小值 fv
% 调用格式：[xm,fv] = LnCPSO(fitness,N,c1,c2,wmax,wmin,M,D) 

function [xm,fv] = LnCPSO(fitness,N,cmax,cmin,w,M,D)

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
    c = cmax - (cmax-cmin)*t/M;
    for i = 1:N
        v(i,:) = w*v(i,:)+c*rand*(y(i,:)-x(i,:))+c*rand*(pg-x(i,:));
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
        