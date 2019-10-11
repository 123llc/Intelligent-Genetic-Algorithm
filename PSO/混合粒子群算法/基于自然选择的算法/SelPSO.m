%% ������Ȼѡ����㷨
% N ������Ŀ
% c1 ѧϰ����1
% c2 ѧϰ����2
% w ����Ȩ��
% M ����������
% D �Ա����ĸ��� 
% xm Ŀ�꺯��ȡ��Сֵ���Ա���ֵ 
% fv Ŀ�꺯������Сֵ 
% fitness ���Ż���Ŀ�꺯�� 
% ���ø�ʽ��[xm,fv] = SelPSO(fitness,N,c1,c2,w,M,D)

function [xm,fv] = SelPSO(fitness,N,c1,c2,w,M,D)

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

for t = 1:M  % M ����������
    for i = 1:N % N ������Ŀ
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
   [sortf,sortx] = sort(fx); %�����е����Ӱ���Ӧֵ����
   exIndex = round((N-1)/2);
   % ����õ�һ�����ӵ�λ���滻������һ��
   x(sortx((N-exIndex+1):N)) = x(sortx(1:exIndex));
   % ����õ�һ�����ӵ��ٶ��滻������һ��
   y(sortx((N-exIndex+1):N)) = v(sortx(1:exIndex));
end
xm = pg';
fv = fitness(pg);     