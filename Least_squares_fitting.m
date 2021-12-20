
%% 
% ��������ֵ����2020�����о������ϻ�ʵ��
% �������� | ��С�������/Least squares fitting
% ʱ��: 2020/10/29
% ѧ��: 20S****** ***
% ----------------------------------------------------------
% 1����������С���˷�����ʵ�����ݡ�
%%
% ������ϵĺ���ģ�ͼ�����ʽ����ߵ����������ݵ�
syms f x;
r = 1;
X = [3 4 5 6 7 8 9];
Y = [2.01 2.98 3.50 5.02 5.47 6.02 7.05];
% ���ú����������
[answer] = Lsf(r,X,Y);
fprintf("����ʽ��С������� ϵ���ֱ�Ϊ\n");
for i = 0:r
    fprintf("  %f  ",answer(i+1));
end
fprintf("\n");
% ����ͼ��Ƚ���ϵĽ��
scatter(X,Y,15,'filled');
hold on;
grid on;
num = 0;
fit_x = zeros(1,ceil((max(X)-min(X))/0.01)+1);
fit_y = zeros(1,ceil((max(X)-min(X))/0.01)+1);
for m = min(X):0.01:max(X)
    num = num + 1;
    fit_x(num) = m;
    for j = 0:r
        fit_y(num) = fit_y(num) + answer(j+1)*fit_x(num)^j;
    end
end
plot(fit_x,fit_y);
%%
% ���庯��������ʽ��С�������ϵ��
function [answer] = Lsf(r,X,Y)
% ���巨����ϵ�����󼰷����Ҳ�����
A = zeros(r+1);
b = zeros(r+1,1);
% ����С����ԭ��/���ƽ���ƽ���������̸���ϵ��/��Ԫ������ֵ
for i = 0:r
    for j = 0:r
        for n = 0:length(X)-1
            A(i+1,j+1) = A(i+1,j+1) + 1*X(n+1)^(i+j);
        end
    end
end
% ����������Ҳ�����
for i = 0:r
    for n = 0:length(X)-1
        b(i+1) = b(i+1) + 1*X(n+1)^i*Y(n+1);
    end
end
% ���ϵ������������ϸ�ʵ���ж���ĸ�˹����Ԫ��ȥ������������Է�����
[answer] = gauss_principal_element_elimination(A,b);
end

%%
%���庯��ʵ�ָ�˹����Ԫ��ȥ��
function [answer] = gauss_principal_element_elimination(A,b)
% ��ȡϵ������Ľ���
[count] = size(A,1);
% ��Ԫ����
for i = 1:count-1
    % ÿһ����Ԫǰ������ѡ��Ԫ
    cursor = i; max = abs(A(i,i));
    for m = i:count
       if abs(A(m,i))>max
           max = abs(A(m,i));
           cursor = m;
       end
    end
    % ����Ԫ�н���
    if cursor ~= i
        row_temp = A(i,:);
        A(i,:) = A(cursor,:);
        A(cursor,:) = row_temp;
        b_temp = b(i);
        b(i) = b(cursor);
        b(cursor) = b_temp;
    end
    % ��Ԫ
    for j = i:count-1
        for n = 1:count-i
            A(i+n,j+1) =  A(i+n,j+1)-A(i,j+1)*A(i+n,i)/A(i,i);
        end
    end
    for n = 1:count-i
        b(i+n) = b(i+n)-b(i)*A(i+n,i)/A(i,i);
    end
end
% �ش�������
answer = zeros(count,1);
for i = 1:count
    answer(count-i+1) = b(count-i+1);
    if i>1
        for j = 1:i-1
            answer(count-i+1) = answer(count-i+1) - answer(count-i+1+j)*A(count-i+1,count-i+1+j);
        end
    end
    answer(count-i+1) = answer(count-i+1)/A(count-i+1,count-i+1);
end
end
%% ------------------END OF THE FILE------------------