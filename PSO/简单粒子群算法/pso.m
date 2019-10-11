%% ��������Ⱥ�㷨
% ������Ŀ��N
% ѧϰ���� c1,c2
% ����Ȩ�� w
% ���������� M 
% �Ա����ĸ��� D
% Ŀ�꺯��ȡֵ��Сʱ���Ա���ֵ xm
% Ŀ�꺯������Сֵ fv

% ���ø�ʽ��[xm,fv] = pso(fitness,N,c1,c2,w,M,D) 

function [xm,fv] = pso(fitness,N,c1,c2,w,M,D)

format long;
% ��ʼ���ٶȺ�λ��
for i = 1:N
    for j = 1:D
        x(i,j) = randn;
        v(i,j) = randn;
    end
end
% p(i) �ǲ�������Ӧ�ȵ�ֵ
% y(i,:)��x(i,:) �����λ�õ�ֵ
for i = 1:N
    p(i) = fitness(x(i,:)); % ��Ӧ��ֵ
    y(i,:) = x(i,:);
end

pg = x(N,:);  % pgΪȫ������
for i = 1:(N-1)
    if fitness(x(i,:) < fitness(pg))
        pg = x(i,:);
    end
end

% �ٶȡ�λ�ø���
for t = 1:M
    for i = 1:N
        v(i,:) = w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        % ��Ӧ��ֵ�ıȽ�
        if fitness(x(i,:)) < p(i)
            y(i,:)=x(i,:);
        end
        % λ����Ϣ�ıȽ�
        if p(i) < fitness(pg)
            pg = y(i,:);
        end
    end
    % ����λ��
    Pbest(t) = fitness(pg);
end
xm = pg';
fv = fitness(pg);
        
