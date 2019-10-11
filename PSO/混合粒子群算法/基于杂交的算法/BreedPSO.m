%% �����ӽ����㷨
% N ������Ŀ
% c1 ѧϰ����1
% c2 ѧϰ����2
% w ����Ȩ��
% Pc �ӽ�����
% Sp �ӽ��صĴ�С����
% M ����������
% D �Ա����ĸ��� 
% xm Ŀ�꺯��ȡ��Сֵ���Ա���ֵ 
% fv Ŀ�꺯������Сֵ 
% fitness ���Ż���Ŀ�꺯�� 
% ���ø�ʽ��[xm,fv] = BreedPSO(fitness,N,c1,c2,w,M,D)

function [xm,fv] = BreedPSO(fitness,N,c1,c2,w,Pc,Sp,M,D)

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
        if fitness(x(i,:)) < p(i)
            p(i) = fitness(x(i,:));
            y(i,:) = x(i,:);
        end
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
        r1 = rand();
        if r1 < Pc % �ӽ�����
            numPool = round(Sp*N); % �ӽ��ش�С
            PoolX = x(1:numPool,:); % �ӽ��������ӵ�λ��
            PoolVX = x(1:numPool,:); % �ӽ��������ӵ��ٶ�
            for i=1:numPool
                seed1 = floor(rand()*(numPool-1)) + 1;
                seed2 = floor(rand()*(numPool-1)) + 1;
                pb = rand();
                % �Ӵ���λ�ü���
                childx1(i,:) = pb*PoolX(seed1,:) + (1-pb)*PoolX(seed2,:);
                % �Ӵ����ٶȼ���
                childv1(i,:) = (PoolVX(seed1,:) +PoolVX(seed2,:))*norm(PoolVX(seed1,:))/...
                    norm(PoolVX(seed1,:) + PoolVX(seed2,:));
            end
            x(1:numPool,:) = childx1; % �Ӵ���λ���滻������λ��
            v(1:numPool,:) = childv1; % �Ӵ����ٶ��滻�������ٶ�
        end
    end
end

xm = pg';
fv = fitness(pg);     