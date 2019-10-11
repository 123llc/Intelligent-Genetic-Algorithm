%% ��ѧϰ���ӵ�����Ⱥ�㷨-ͬ���仯��ѧϰ����
% ������Ŀ��N
% ѧϰ���ӵ����ֵ��cmax
% ѧϰ���ӵ���Сֵ��cmin
% ���Ȩ�صķ���
% ���������� M 
% �Ա����ĸ��� D
% Ŀ�꺯��ȡֵ��Сʱ���Ա���ֵ xm
% Ŀ�꺯������Сֵ fv
% ���ø�ʽ��[xm,fv] = LnCPSO(fitness,N,c1,c2,wmax,wmin,M,D) 

function [xm,fv] = LnCPSO(fitness,N,cmax,cmin,w,M,D)

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
    c = cmax - (cmax-cmin)*t/M;
    for i = 1:N
        v(i,:) = w*v(i,:)+c*rand*(y(i,:)-x(i,:))+c*rand*(pg-x(i,:));
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
        