%% 基于杂交的算法
% N 粒子数目
% c1 学习因子1
% c2 学习因子2
% lamda 退火常数
% M 最大迭代次数
% D 自变量的个数 
% xm 目标函数取最小值的自变量值 
% fv 目标函数的最小值 
% fitness 待优化的目标函数 
% 调用格式：[xm,fv] = SimuAPSO(fitness,N,c1,c2,lamda,M,D)

function [xm,fv] = SimuAPSO(fitness,N,c1,c2,lamda,M,D)

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

T = fitness(pg)/log(5) % 初始化温度
for t = 1:M  % M 最大迭代次数
    groupFit = fitness(pg);
    for i = 1:N % N 粒子数目
        Tfit(i) = exp(-(p(i) - groupFit)/T);
    end
    SumTfit = sum(Tfit);
    Tfit = Tfit/SumTfit;
    pBet = rand();
    for i=1:N % 当前温度下各个pi的适应值   N 粒子数目
        ComFit(i) = sum(Tfit(1:i));
        if pBet <= ComFit(i)
            pg_plus = x(i);
            break;
        end
    end
    C = c1 + c2;
    ksi = 2/abs(2-C-sqrt(C^2-4*C)); % 速度压缩因子
    for i=1:N
        v(i,:)=ksi*(v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg_plus-x(i,:)));
        x(i,:) = x(i,:) + v(i,:);
        if fitness(x(i,:)) < p(i)
            p(i)=fitness(x(i,:));
            y(i,:) = x(i,:);
        end
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
    T = T*lamda; % 退温操作
end     
xm = pg';
fv = fitness(pg);     