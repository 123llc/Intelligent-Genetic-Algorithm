%% 二阶振荡粒子群算法
% N 粒子数目
% c1 学习因子1
% c2 学习因子2
% w 惯性权重
% xmax 自变量搜索域的最大值
% xmin 自变量搜索域的最小值
% M 最大迭代次数
% MaxC 混沌搜索的最大步数
% D 自变量的个数 
% xm 目标函数取最小值的自变量值 
% fv 目标函数的最小值 
% fitness 待优化的目标函数 
% 调用格式：[xm,fv] = CLSPSO(fitness,N,w,c1,c2,xmax,xmin,M,MaxC,D)

function [xm,fv] = CLSPSO(fitness,N,w,c1,c2,xmax,xmin,M,MaxC,D)

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
        fv(i) = fitness(x(i,:));
    end
    [sort_fv,index] = sort(fv); %sort_fv从小到大（左到右），index以前的大小的排序值
    Nbest = floor(N*0.2); %保留群体中20%的最佳粒子
    for n=1:Nbest % 对群体中20%的最佳粒子进行混沌搜索
        tmpx = x(index(n),:);
        for k=1:MaxC % 混沌搜索的最大步长
            for dim = 1:D % 混沌搜索的迭代公式
                cx(dim) = (tmpx(1,dim) - xmin(dim))/(tmpx(1,dim) - xmax(dim));
                cx(dim) = 4*cx(dim)*(1-cx(dim));
                tmpx(1,dim) = tmpx(1,dim) + cx(dim)*(xmax(dim) - xmin(dim));
            end
            fcs = fitness(tmpx);
            if fcs < sort_fv(n) % 对混沌搜索后的决策变量值进行评估
                x(index(n),:) = tmpx;
                break;
            end
        end
            x(index(n),:) = tmpx;
    end
    r = rand();
    for s=1:D % 收缩搜索区域
        xmin(s) = max(xmin(s), pg(s)-r*(xmax(s)-xmin(s)));
        xmax(s) = min(xmax(s), pg(s)+r*(xmax(s)-xmin(s)));
    end
    x(1:Nbest,:) = x(index(1:Nbest),:);
    for i=(Nbest+1):N % 随机产生剩余的80%微粒
        for j=1:D
            x(i,j) = xmin(j) + rand*(xmax(j)-xmin(j)); % 随机初始化位置
            v(i,j) = randn; % 随机初始化速度
        end
    end
    Pbest(t) = fitness(pg);
    for i=1:N
        if fitness(x(i,:)) < p(i)
            p(i) = fitness(x(i,:));
            y(i,:) = x(i,:);
        end
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
end
xm = pg';
fv = fitness(pg);     