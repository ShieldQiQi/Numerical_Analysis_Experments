
%% 
% 哈工大数值分析2020年秋研究生，上机实验
% 第三部分 | 最小二乘拟合/Least squares fitting
% 时间: 2020/10/29
% 学生: 20S****** ***
% ----------------------------------------------------------
% 1、【利用最小二乘法处理实验数据】
%%
% 定义拟合的函数模型及多项式的最高的幂数和数据点
syms f x;
r = 1;
X = [3 4 5 6 7 8 9];
Y = [2.01 2.98 3.50 5.02 5.47 6.02 7.05];
% 调用函数进行拟合
[answer] = Lsf(r,X,Y);
fprintf("多项式最小二乘拟合 系数分别为\n");
for i = 0:r
    fprintf("  %f  ",answer(i+1));
end
fprintf("\n");
% 绘制图像比较拟合的结果
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
% 定义函数求解多项式最小二乘拟合系数
function [answer] = Lsf(r,X,Y)
% 定义法方程系数矩阵及方程右侧向量
A = zeros(r+1);
b = zeros(r+1,1);
% 由最小二乘原则/最佳平方逼近求出法方程各个系数/多元函数求极值
for i = 0:r
    for j = 0:r
        for n = 0:length(X)-1
            A(i+1,j+1) = A(i+1,j+1) + 1*X(n+1)^(i+j);
        end
    end
end
% 算出法方程右侧向量
for i = 0:r
    for n = 0:length(X)-1
        b(i+1) = b(i+1) + 1*X(n+1)^i*Y(n+1);
    end
end
% 求出系数矩阵后，利用上个实验中定义的高斯列主元消去法求解非齐次线性方程组
[answer] = gauss_principal_element_elimination(A,b);
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