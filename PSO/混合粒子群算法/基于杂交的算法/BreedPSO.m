%% 基于杂交的算法
% N 粒子数目
% c1 学习因子1
% c2 学习因子2
% w 惯性权重
% Pc 杂交概率
% Sp 杂交池的大小比例
% M 最大迭代次数
% D 自变量的个数 
% xm 目标函数取最小值的自变量值 
% fv 目标函数的最小值 
% fitness 待优化的目标函数 
% 调用格式：[xm,fv] = BreedPSO(fitness,N,c1,c2,w,M,D)

function [xm,fv] = BreedPSO(fitness,N,c1,c2,w,Pc,Sp,M,D)

format long; % 设置输出格式
for i=1:N
    for j=1:D
        x(i,j) = randn; % 随机初始化位置
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

for t = 1:M  % M 最大迭代次数
    for i = 1:N % N 粒子数目
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        if fitness(x(i,:)) < p(i)
            p(i) = fitness(x(i,:));
            y(i,:) = x(i,:);
        end
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
        r1 = rand();
        if r1 < Pc % 杂交概率
            numPool = round(Sp*N); % 杂交池大小
            PoolX = x(1:numPool,:); % 杂交池中粒子的位置
            PoolVX = x(1:numPool,:); % 杂交池中粒子的速度
            for i=1:numPool
                seed1 = floor(rand()*(numPool-1)) + 1;
                seed2 = floor(rand()*(numPool-1)) + 1;
                pb = rand();
                % 子代的位置计算
                childx1(i,:) = pb*PoolX(seed1,:) + (1-pb)*PoolX(seed2,:);
                % 子代的速度计算
                childv1(i,:) = (PoolVX(seed1,:) +PoolVX(seed2,:))*norm(PoolVX(seed1,:))/...
                    norm(PoolVX(seed1,:) + PoolVX(seed2,:));
            end
            x(1:numPool,:) = childx1; % 子代的位置替换父代的位置
            v(1:numPool,:) = childv1; % 子代的速度替换父代的速度
        end
    end
end

xm = pg';
fv = fitness(pg);     