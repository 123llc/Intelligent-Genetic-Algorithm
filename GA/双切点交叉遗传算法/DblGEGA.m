%% 双切点遗传算法
% a 自变量下届
% b 自变量上届
% NP 种群大小
% NG 最大进化代数
% Pc 杂交概率
% Pm 变异概率
% eps 自变量离散精度
% xm 目标函数取最小值时的自变量值
% fv 目标函数的最小值
function [xv,fv] = DblGEGA(fitness,a,b,NP,NG,Pc,Pm,eps)

L = ceil(log2((b-a)/eps+1)); % 根据离散精度，确定二进制需要的码长
x = zeros(NP,L); % NP种群大小，L编码长度，生成一个矩阵

% 初始化
for i = 1:NP
    x(i,:) = Initial(L); % 种群初始化
    fx(i) = fitness(Dec(a,b,x(i,:),L)); % 个体适应值
end

for k=1:NG
    sumfx = sum(fx); %所有个体适应度之和
    Px = fx/sumfx;
    PPx = 0;
    PPx(1) = Px(1);
    for i=2:NP
        PPx(i) = PPx(i-1) + Px(i);
    end
    for i=1:NP
        sita = rand();
        for n = 1:NP
            if sita <= PPx(n)
                SelFather = n; % 根据轮盘赌策略确定的父亲
                break;
            end
        end
        Selmother = floor(rand()*(NP-1))+1; % 随机选择母亲
        posCut1 = round(rand()*(L-2))+1;    % 随机确定交叉位置 1
        posCut2 = round(rand()*(L-2))+1;    % 随机确定交叉位置 2
        while posCut2 == posCut1 % 保证两个交叉位置不同
            posCut2 = round(rand()*(L-2))+1;
        end
        posCutmax = max(posCut1,posCut2); % 确定两个交叉位置中的较大值
        posCutmin = min(posCut1,posCut2); % 确定两个交叉位置中的较小值
        r1 = rand();
        if r1 <= Pc % 交叉
            nx(i,1:posCutmin) = x(SelFather,1:posCutmin); % 交换posCutmin到posCutmax之间的部分
            nx(i,(posCutmin+1):posCutmax) = x(Selmother,(posCutmin+1):posCutmax);
            if posCutmax<L
                nx(i,(posCutmax+1):L) = x(SelFather,(posCutmax+1):L);
            end
            r2 = rand();
            if r2 <= Pm % 变异
                posMut = round(rand()*(L-1) + 1);
                nx(i,posMut) = ~nx(i,posMut);
            end
        else
            nx(i,:) = x(SelFather,:);
        end
    end
    x = nx;
    for i=1:NP
        fx(i) = fitness(Dec(a,b,x(i,:),L));
    end
end

fv = -inf; % matlab输出的结果很大时的一种声明
for i=1:NP
    fitx = fitness(Dec(a,b,x(i,:),L)); % 取个体中的最好值作为最终结果
    if fitx > fv
        fv = fitx;
        xv = Dec(a,b,x(i,:),L);
    end
end
    
function result = Initial(length) % 初始化函数
for i=1:length
    r = rand();
    result(i) = round(r);
end

function y = Dec(a,b,x,L) % 二进制编码转为十进制编码
base = 2.^((L-1):-1:0);
y = dot(base,x);
y = a+y*(b-a)/(2^L-1); 