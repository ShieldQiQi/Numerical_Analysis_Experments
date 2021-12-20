
%% 
% ��������ֵ����2020�����о������ϻ�ʵ��
% �ڶ����� | ���Է��������/��˹����Ԫ��ȥ��
% ʱ��: 2020/10/22
% ѧ��: 20S****** ***
% ----------------------------------------------------------
% 1������˹��ȥ����
% 2������˹����Ԫ��ȥ����

% ����˹��ȥ�������Է�����ʱ�������ԪԪ�ص���0������ȥ���޷�������������ԪԪ�ؽӽ���0������ʹ����ȥ�������²��ȶ�����
% ��ʱ��Ҫʹ�ø�˹����Ԫ��ȥ��
%%
% ����Ҫ���ķ�����һ
A = [10^(-8) 2 3; -1 3.712 4.623; -2 1.072 5.643];
b = [1 2 3];
% ����Ҫ���ķ������
C = [4 -2 4; -2 17 10; -4 10 9];
d = [10 3 7];
%%
% �ֱ��ø�˹����Ԫ��ȥ���͸�˹��ȥ����ⷽ����
[answer] = gauss_elimination(A,b);
fprintf("��˹�� ����һ answer = [ %10f %f %f ]\n",answer(1),answer(2),answer(3));
[answer] = gauss_principal_element_elimination(A,b);
fprintf("��˹����Ԫ�� ����һ answer = [ %10f %f %f ]\n",answer(1),answer(2),answer(3));
[answer] = gauss_elimination(C,d);
fprintf("��˹�� ���̶� answer = [ %f %f %f ]\n",answer(1),answer(2),answer(3));
[answer] = gauss_principal_element_elimination(C,d);
fprintf("��˹����Ԫ�� ���̶� answer = [ %f %f %f ]\n",answer(1),answer(2),answer(3));
%%
%���庯��ʵ�ָ�˹��ȥ��
function [answer] = gauss_elimination(A,b)
% ��ȡϵ������Ľ���
[count] = size(A,1);
% ��Ԫ����
for i = 1:count-1
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