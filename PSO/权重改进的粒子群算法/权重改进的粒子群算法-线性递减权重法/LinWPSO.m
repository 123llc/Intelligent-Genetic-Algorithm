%% Ȩ�ظĽ�������Ⱥ�㷨
% ������Ŀ��N
% ѧϰ���� c1,c2
% ����Ȩ�� w
% ���Ȩ�� wmax
% ��СȨ�� wmin
% ���������� M 
% �Ա����ĸ��� D
% Ŀ�꺯��ȡֵ��Сʱ���Ա���ֵ xm
% Ŀ�꺯������Сֵ fv

% ���ø�ʽ��[xm,fv] = LinYSPSO(fitness,N,c1,c2,wmax,wmin,M,D) 

function [xm,fv] = YSPSO(fitness,N,c1,c2,wmax,wmin,M,D)

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
    for i = 1:N
        w = wmax - (t-1)*(wmax-wmin)/(M-1); % Ȩ�����Եݼ�
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        if fitness(x(i,:)) < p(i)
            p(i) =fitness(x(i,:));
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
        