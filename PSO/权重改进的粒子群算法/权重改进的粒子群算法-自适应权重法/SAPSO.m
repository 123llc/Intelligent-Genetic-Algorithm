%% 自适应权重法
% 粒子数目：N
% 学习因子 c1,c2
% 最大权重 wmax
% 最小权重 wmin
% 最大迭代次数 M 
% 自变量的个数 D
% 目标函数取值最小时的自变量值 xm
% 目标函数的最小值 fv

% 调用格式：[xm,fv] = SAPSO(fitness,N,c1,c2,wmax,wmin,M,D) 

function [xm,fv] = SAPSO(fitness,N,c1,c2,wmax,wmin,M,D)

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
    for j = 1:N
        fv(j) = fitness(x(j,:));
    end
    fvag = sum(fv)/N; % 适应度平均值
    fmin = min(fv);   % 适应度最小值
    for i = 1:40    
        if fv(i) <= fvag % 自适应权重计算公式
            w = wmin + (fv(i) - fmin)*(wmax - wmin)/(fvag - fmin);
        else
            w = wmax;
        end
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
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
        