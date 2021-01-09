
%% 
% 哈工大数值分析2020年秋研究生，上机实验
% 第二部分 | 线性方程组求解/高斯列主元消去法
% 时间: 2020/10/22
% 学生: 20S****** ***
% ----------------------------------------------------------
% 1、【高斯消去法】
% 2、【高斯列主元消去法】

% 若高斯消去法解线性方程组时，如果主元元素等于0，则消去法无法继续，或者主元元素接近于0，继续使用消去法将导致不稳定现象，
% 此时需要使用高斯列主元消去法
%%
% 定义要求解的方程组一
A = [10^(-8) 2 3; -1 3.712 4.623; -2 1.072 5.643];
b = [1 2 3];
% 定义要求解的方程组二
C = [4 -2 4; -2 17 10; -4 10 9];
d = [10 3 7];
%%
% 分别用高斯列主元消去法和高斯消去法求解方程组
[answer] = gauss_elimination(A,b);
fprintf("高斯法 方程一 answer = [ %10f %f %f ]\n",answer(1),answer(2),answer(3));
[answer] = gauss_principal_element_elimination(A,b);
fprintf("高斯列主元法 方程一 answer = [ %10f %f %f ]\n",answer(1),answer(2),answer(3));
[answer] = gauss_elimination(C,d);
fprintf("高斯法 方程二 answer = [ %f %f %f ]\n",answer(1),answer(2),answer(3));
[answer] = gauss_principal_element_elimination(C,d);
fprintf("高斯列主元法 方程二 answer = [ %f %f %f ]\n",answer(1),answer(2),answer(3));
%%
%定义函数实现高斯消去法
function [answer] = gauss_elimination(A,b)
% 获取系数矩阵的阶数
[count] = size(A,1);
% 消元过程
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
% 回代求解过程
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
%定义函数实现高斯列主元消去法
function [answer] = gauss_principal_element_elimination(A,b)
% 获取系数矩阵的阶数
[count] = size(A,1);
% 消元过程
for i = 1:count-1
    % 每一次消元前进行列选主元
    cursor = i; max = abs(A(i,i));
    for m = i:count
       if abs(A(m,i))>max
           max = abs(A(m,i));
           cursor = m;
       end
    end
    % 列主元行交换
    if cursor ~= i
        row_temp = A(i,:);
        A(i,:) = A(cursor,:);
        A(cursor,:) = row_temp;
        b_temp = b(i);
        b(i) = b(cursor);
        b(cursor) = b_temp;
    end
    % 消元
    for j = i:count-1
        for n = 1:count-i
            A(i+n,j+1) =  A(i+n,j+1)-A(i,j+1)*A(i+n,i)/A(i,i);
        end
    end
    for n = 1:count-i
        b(i+n) = b(i+n)-b(i)*A(i+n,i)/A(i,i);
    end
end
% 回代求解过程
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