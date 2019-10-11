%% ����������Ⱥ�㷨
% ������Ŀ��N
% ѧϰ����1��c1
% ѧϰ����2��c2
% ����Ȩ�� w
% ���������� M 
% �Ա����ĸ��� D
% Ŀ�꺯��ȡ��Сֵ���Ա���ֵ xm
% Ŀ�꺯������Сֵ fv
% ���Ż���Ŀ�꺯�� fitness
% ���ø�ʽ��[xm,fv] = AsyLnCPSO(fitness,N,c1max,clmin,c2max,c2min,w,M,D) 

function [xm,fv] = SecVibratPSO(fitness,N,w,c1,c2,M,D)

format long; % ���������ʽ
for i=1:N
    for j=1:D
        x(i,j) = randn; % �����ʼ��λ��
        x1(i,j) = randn; % �����ʼ��λ�ã����ڱ��������ϴε�λ��
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

for t = 1:M
    for i = 1:N
        phi1 = c1*rand();
        phi2 = c2*rand();
        if t < M/2
            ks1 = (2*sqrt(phi1)-1)*rand()/phi1;
            ks2 = (2*sqrt(phi2)-1)*rand()/phi2;
        else
            ks1 = (2*sqrt(phi1)-1)*(1+rand())/phi1;
            ks2 = (2*sqrt(phi2)-1)*(1+rand())/phi2;
        end
        v(i,:) = w*v(i,:)+phi1*(y(i,:)-(1+ks1)*x(i,:)+ks1*x1(i,:))+...
            phi2*(pg-(1+ks2)*x(i,:)+ks1*x1(i,:));
        x1(i,:) = x(i,:);
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
        