%% 多变异位自适应遗传算法
% a 自变量下届
% b 自变量上届
% NP 种群大小
% NG 最大进化代数
% Pc1 杂交常数1
% Pc2 杂交常数2
% Pm1 变异常数1
% Pm2 变异常熟2
% eps 自变量离散精度
% xm 目标函数取最小值时的自变量值
% fv 目标函数的最小值
function [xv,fv] = MMAdapGA(fitness,a,b,NP,NG,Pc1,Pc2,Pm1,Pm2,eps)

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
        posCut = round(rand()*(L-2))+1;    % 随机确定交叉位置
        favg = sumfx/NP; % 群体平均适应度
        fmax = max(fx); % 群体最大适应值
        Fitness_f = fx(SelFather); % 交叉的父亲适应值
        Fitness_m = fx(Selmother); % 交叉的母亲适应值
        Fm = max(Fitness_f,Fitness_m); % 交叉双方较大的适应值
        if Fm >= favg % 自适应交叉概率
            Pc = Pc1*(fmax - Fm)/(fmax - favg);
        else
            Pc = Pc2;
        end
        r1 = rand();
        if r1 <= Pc % 交叉
            nx(i,1:posCut) = x(SelFather,1:posCut);
            nx(i,(posCut+1):L) = x(Selmother,(posCut+1):L);
            fmu = fitness(Dec(a,b,nx(i,:),L));
            if fmu>=favg % 自适应变异概率
                Pm = Pm1*(fmax - fmu)/(fmax - favg);
            else
                Pm = Pm2;
            end
            r2 = rand();
            if r2 <= Pm %变异 计算numMut表示有多少个变异位
                if fmax ~= favg
                    numMut = round(L*(fmax-fmu)/(fmax - favg)/3); % 随机产生变异位置的个数
                else
                    numMut = 1; % 如果群体中的最大值等于最小值，只产生一个变异位
                end
                intArr = 1:L;
                posMut = zeros(1,numMut);
                for j=1:numMut % 随机产生变异位置编号
                    tmp = round(rand()*(L-j) + 1);
                    posMut(j) = intArr(tmp);
                    intArr = [intArr(1:(tmp-1)) intArr((tmp+1):(L-j+1))]; % 动态的更新数组，去掉取走的值
                end
                for j=1:numMut % 对变异位进行取反
                    nx(i,posMut(j)) = ~ nx(i,posMut(j));
                end
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
end
fv = -inf;
for i=1:NP
    fitx = fitness(Dec(a,b,x(i,:),L));
    if fitx > fv
        fv = fitx; % 取个体中最好值作为最终结果
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