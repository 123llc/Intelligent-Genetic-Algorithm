%% ˳��ѡ���Ŵ��㷨
% a �Ա����½�
% b �Ա����Ͻ�
% NP ��Ⱥ��С
% NG ����������
% q ��ø����ѡ�����
% Pc �ӽ�����
% Pm �������
% eps �Ա�����ɢ����
% xm Ŀ�꺯��ȡ��Сֵʱ���Ա���ֵ
% fv Ŀ�꺯������Сֵ
function [xv,fv] = SBOGA(fitness,a,b,NP,NG,q,Pc,Pm,eps)

L = ceil(log2((b-a)/eps+1)); % ������ɢ���ȣ�ȷ����������Ҫ���볤
x = zeros(NP,L); % NP��Ⱥ��С��L���볤�ȣ�����һ������\

% ��ʼ��
for i = 1:NP
    x(i,:) = Initial(L); % ��Ⱥ��ʼ��
    fx(i) = fitness(Dec(a,b,x(i,:),L)); % ������Ӧֵ
end

for i=1:NG % NG ����������
    [sortf,sortx] = sort(fx);% ��Ӧֵ����
    x = x(sortx,:);
    fx = fx(sortx);
    for i=1:NP % �̶�ѡ����� 
        Px(i) = (1-q)^(NP-i)*q/(1-(1-q)^NP);
    end
    PPx = 0;
    PPx(1) = Px(1);
    for i=2:NP % �������̶Ĳ��Եĸ����ۼ�
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
        posCut = floor(rand()*(L-2))+1; % ���ȷ�������
        r1 = rand();
        if r1 <= Pc % ����
            nx(i,1:posCut) = x(SelFather,1:posCut);
            nx(i,(posCut+1):L) = x(Selmother,(posCut+1):L);
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
        fx(i) = fitness(Dec(a,b,x(i,:),L)); % �Ӵ���Ӧֵ
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