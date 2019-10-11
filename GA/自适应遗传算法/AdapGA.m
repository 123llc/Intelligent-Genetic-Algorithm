%% ����Ӧ�Ŵ��㷨
% a �Ա����½�
% b �Ա����Ͻ�
% NP ��Ⱥ��С
% NG ����������
% Pc1 �ӽ�����1
% Pc2 �ӽ�����2
% Pm1 ���쳣��1
% Pm2 ���쳣��2
% eps �Ա�����ɢ����
% xm Ŀ�꺯��ȡ��Сֵʱ���Ա���ֵ
% fv Ŀ�꺯������Сֵ
function [xv,fv] = AdapGA(fitness,a,b,NP,NG,Pc1,Pc2,Pm1,Pm2,eps)

L = ceil(log2((b-a)/eps+1)); % ������ɢ���ȣ�ȷ����������Ҫ���볤
x = zeros(NP,L); % NP��Ⱥ��С��L���볤�ȣ�����һ������

% ��ʼ��
for i = 1:NP
    x(i,:) = Initial(L); % ��Ⱥ��ʼ��
    fx(i) = fitness(Dec(a,b,x(i,:),L)); % ������Ӧֵ
end

for k=1:NG
    sumfx = sum(fx); %���и�����Ӧ��֮��
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
                SelFather = n; % �������̶Ĳ���ȷ���ĸ���
                break;
            end
        end
        Selmother = floor(rand()*(NP-1))+1; % ���ѡ��ĸ��
        posCut = floor(rand()*(L-2))+1;     % ���ȷ�������
        favg = sumfx/NP;                    % ��Ⱥƽ����Ӧֵ
        fmax = max(fx);                     % ��Ⱥ�����Ӧֵ
        Fitness_f = fx(Selmother);          % ����ĸ�����Ӧֵ
        Fitness_m = fx(Selmother);          % �����ĸ����Ӧֵ
        Fm = max(Fitness_f,Fitness_m);      % ����Ӧ�������
        if Fm >= favg
            Pc = Pc1*(fmax - Fm)/(fmax - favg);
        else
            Pc = Pc2;
        end
        r1 = rand();
        if r1 <= Pc % ����
            nx(i,1:posCut) = x(SelFather,1:posCut);
            nx(i,(posCut+1):L) = x(Selmother,(posCut+1):L);
            fmu = fitness(Dec(a,b,nx(i,:),L));
            if fmu >= favg % ����Ӧ�������
                Pm = Pm1*(fmax - fmu)/(fmax - favg);
            else
                Pm = Pm2;
            end
            r2 = rand();
            if r2 <= Pm % ����
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

fv = -inf; % matlab����Ľ���ܴ�ʱ��һ������
for i=1:NP
    fitx = fitness(Dec(a,b,x(i,:),L)); % ȡ�����е����ֵ��Ϊ���ս��
    if fitx > fv
        fv = fitx;
        xv = Dec(a,b,x(i,:),L);
    end
end
    
function result = Initial(length) % ��ʼ������
for i=1:length
    r = rand();
    result(i) = round(r);
end

function y = Dec(a,b,x,L) % �����Ʊ���תΪʮ���Ʊ���
base = 2.^((L-1):-1:0);
y = dot(base,x);
y = a+y*(b-a)/(2^L-1); 