%% 二阶振荡粒子群算法
% 粒子数目：N
% 学习因子1：c1
% 学习因子2：c2
% 惯性权重 w
% 最大迭代次数 M 
% 自变量的个数 D
% 目标函数取最小值的自变量值 xm
% 目标函数的最小值 fv
% 待优化的目标函数 fitness
% 调用格式：[xm,fv] = AsyLnCPSO(fitness,N,c1max,clmin,c2max,c2min,w,M,D) 

function [xm,fv] = SecVibratPSO(fitness,N,w,c1,c2,M,D)

format long; % 设置输出格式
for i=1:N
    for j=1:D
        x(i,j) = randn; % 随机初始化位置
        x1(i,j) = randn; % 随机初始化位置，用于保存粒子上次的位置
        v(i,j) = randn; % 随机初始化速度，用于保存粒子当前的位置
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
    for i = 1:N
        phi1 = c1*rand();
        phi2 = c2*rand();
        if t < M/2
            ks1 = (2*sqrt(phi1)-1)*rand()/phi1;
            ks2 = (2*sqrt(phi2)-1)*rand()/phi2;
        else
            ks1 = (2*sqrt(phi1)-1)*(1+rand())/phi1;
            ks2 = (2*sqrt(phi2)-1)*(1+rand())/phi2;
        end
        v(i,:) = w*v(i,:)+phi1*(y(i,:)-(1+ks1)*x(i,:)+ks1*x1(i,:))+...
            phi2*(pg-(1+ks2)*x(i,:)+ks1*x1(i,:));
        x1(i,:) = x(i,:);
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
        