%% 基于自然选择的算法
% N 粒子数目
% c1 学习因子1
% c2 学习因子2
% w 惯性权重
% M 最大迭代次数
% D 自变量的个数 
% xm 目标函数取最小值的自变量值 
% fv 目标函数的最小值 
% fitness 待优化的目标函数 
% 调用格式：[xm,fv] = SelPSO(fitness,N,c1,c2,w,M,D)

function [xm,fv] = SelPSO(fitness,N,c1,c2,w,M,D)

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
        fx(i) = fitness(x(i,:));
        if fx(i) < p(i)
            p(i) = fx(i);
            y(i,:) = x(i,:);
        end
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
   [sortf,sortx] = sort(fx); %将所有的粒子按适应值排序
   exIndex = round((N-1)/2);
   % 将最好的一半粒子的位置替换掉最差的一半
   x(sortx((N-exIndex+1):N)) = x(sortx(1:exIndex));
   % 将最好的一半粒子的速度替换掉最差的一半
   y(sortx((N-exIndex+1):N)) = v(sortx(1:exIndex));
end
xm = pg';
fv = fitness(pg);     