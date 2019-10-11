%% ��ѧϰ���ӵ�����Ⱥ�㷨-�첽�仯��ѧϰ����
% ������Ŀ��N
% ѧϰ����1�����ֵ��c1max
% ѧϰ����1����Сֵ��c1min
% ѧϰ����2�����ֵ��c2max
% ѧϰ����2����Сֵ��c2min
% ����Ȩ�� w
% ���������� M 
% �Ա����ĸ��� D
% Ŀ�꺯��ȡֵ��Сʱ���Ա���ֵ xm
% Ŀ�꺯������Сֵ fv
% ���ø�ʽ��[xm,fv] = AsyLnCPSO(fitness,N,c1max,clmin,c2max,c2min,w,M,D) 

function [xm,fv] = AsyLnCPSO(fitness,N,c1max,c1min,c2max,c2min,w,M,D)

format long; % ���������ʽ
for i=1:N
    for j=1:D
        x(i,j) = randn; %��ʼ���ٶȺ�λ��
        v(i,j) = randn;
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

for t = 1:M
    c1 = c1max - (c1max - c1min)*t/M; %���Ա仯��ѧϰ����1
    c2 = c2max - (c2max - c2min)*t/M; %���Ա仯��ѧУ����2
    for i = 1:N
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        if fitness(x(i,:)) < p(i)
            p(i) = fitness(x(i,:));
            y(i,:) = x(i,:);
        end
        % ÿ��ѭ����Ҫ�Ƚ�
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
end
% '����ת�õ���˼
xm = pg';
fv = fitness(pg);
        