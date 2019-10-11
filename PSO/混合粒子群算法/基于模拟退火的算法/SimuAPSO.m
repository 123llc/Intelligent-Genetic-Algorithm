%% �����ӽ����㷨
% N ������Ŀ
% c1 ѧϰ����1
% c2 ѧϰ����2
% lamda �˻���
% M ����������
% D �Ա����ĸ��� 
% xm Ŀ�꺯��ȡ��Сֵ���Ա���ֵ 
% fv Ŀ�꺯������Сֵ 
% fitness ���Ż���Ŀ�꺯�� 
% ���ø�ʽ��[xm,fv] = SimuAPSO(fitness,N,c1,c2,lamda,M,D)

function [xm,fv] = SimuAPSO(fitness,N,c1,c2,lamda,M,D)

format long; % ���������ʽ
for i=1:N
    for j=1:D
        x(i,j) = randn; % �����ʼ��λ��
        v(i,j) = randn; % �����ʼ���ٶȣ����ڱ������ӵ�ǰ��λ��
    end
end

for i=1:N
    p(i) = fitness(x(i,:)); % p(i) �ǲ�������Ӧ�ȵ�ֵ
    y(i,:) = x(i,:); % y(i,:)��x(i,:) �����λ�õ�ֵ
end

pg = x(N,:);  % pgΪȫ������
for i=1:(N-1)
    if fitness(x(i,:))<fitness(pg)
        pg = x(i,:);
    end
end

T = fitness(pg)/log(5) % ��ʼ���¶�
for t = 1:M  % M ����������
    groupFit = fitness(pg);
    for i = 1:N % N ������Ŀ
        Tfit(i) = exp(-(p(i) - groupFit)/T);
    end
    SumTfit = sum(Tfit);
    Tfit = Tfit/SumTfit;
    pBet = rand();
    for i=1:N % ��ǰ�¶��¸���pi����Ӧֵ   N ������Ŀ
        ComFit(i) = sum(Tfit(1:i));
        if pBet <= ComFit(i)
            pg_plus = x(i);
            break;
        end
    end
    C = c1 + c2;
    ksi = 2/abs(2-C-sqrt(C^2-4*C)); % �ٶ�ѹ������
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
    T = T*lamda; % ���²���
end     
xm = pg';
fv = fitness(pg);     