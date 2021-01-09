
%% 
% 哈工大数值分析2020年秋研究生，上机实验
% 第一部分 | 非线性方程组求解
% 时间: 2020/10/15
% 学生: 20S****** ***
% ----------------------------------------------------------
% 1、【二分法】
% 2、【牛顿法】
% 3、【割线法】
% 4、【改进的牛顿法】
% 5、【拟牛顿法】
%%
% define the algorithm for which be used to solve the nolinear equation
% where the variable can be 0、1、2、3、4 , each one corresponds to the algorithm above

algorithm_index = 4;

%%
% 方法一：二分法
% 题目：用二分法计算方程 sin(x)-pow(x,2)/2 = 0 在(1,2)内的根的近似值，要求ε=0.5*10^(-5)
    
if algorithm_index == 0
    % 易知函数equation(x)在区间(1,2)内连续可导，且方程的根在区间(1,2)内存在且唯一
    a = 1;
    b = 2;
    count = 0;
    % 定义允许的误差
    delta = 0.5*10^(-5);
    % 定义函数原型，方便用于不同的方程，提高程序应用的普遍性
    syms f x;
    f = sin(x)-x^2/2;
    % 求解迭代解以及迭代次数
    [answer,count] = dichotomy(a,b,count,delta,f);
    fprintf("二分法求非线性方程，解为 %f 迭代了 %d 次\n",answer,count);
end
%%
% 方法二：牛顿法
% 题目：用牛顿法求解下列非线性方程的根，题目见实验报告册

if algorithm_index == 1
    % 定义3个方程的初始值
    x1_initial = 0.5;
    x2_initial = 1;
    x3_initial_1 = 0.45;
    x3_initial_2 = 0.65;
    % 定义允许的误差以及最大的迭代次数
    delta = 0.5*10^(-5);
    N = 100;
    % 定义函数原型，方便用于不同的方程，提高程序应用的普遍性
    syms f x;
    % 求解方程1
    count = 0;
    f = x*exp(x)-1;
    [answer,count] = newton(x1_initial,count,delta,f,N);
    fprintf("牛顿法求非线性方程1，初始值为 %f 解为 %f 迭代了 %d 次\n",x1_initial,answer,count);
    % 求解方程2
    count = 0;
    f = x^3-x-1;
    [answer,count] = newton(x2_initial,count,delta,f,N);
    fprintf("牛顿法求非线性方程2，初始值为 %f 解为 %f 迭代了 %d 次\n",x2_initial,answer,count);
    % 求解方程3
    count = 0;
    f = (x-1)^2*(2*x-1);
    [answer,count] = newton(x3_initial_1,count,delta,f,N);
    fprintf("牛顿法求非线性方程3，初始值为 %f 解为 %f 迭代了 %d 次\n",x3_initial_1,answer,count);
    count = 0;
    [answer,count] = newton(x3_initial_2,count,delta,f,N);
    fprintf("牛顿法求非线性方程3，初始值为 %f 解为 %f 迭代了 %d 次\n",x3_initial_2,answer,count);
end
%%
% 方法三：割线法/多点迭代法
% 题目：用割线法求解下列非线性方程的根，题目见实验报告册

if algorithm_index == 2
    % 定义迭代初始点
    x0_initial = 0.4; 
    x1_initial = 0.6; 
    % 定义允许的误差以及最大的迭代次数
    delta = 0.5*10^(-5);
    N = 100;
    count = 0;
    % 定义函数原型，方便用于不同的方程，提高程序应用的普遍性
    syms f x;
    f = x*exp(x)-1;
    % 求解方程  
    [answer,count] = secant(x0_initial,x1_initial,count,delta,f,N);
    fprintf("割线法求非线性方程，初始值x0为 %f 初始值x1为 %f 解为 %f 迭代了 %d 次\n",x0_initial,x1_initial,answer,count);
end
%%
% 方法四：改进的牛顿法
% 题目：用改进的牛顿法求解下列非线性方程的根，题目见实验报告册

if algorithm_index == 3
    % 定义迭代初始点
    x_initial = 0.55; 
    % 定义允许的误差以及最大的迭代次数
    delta = 0.5*10^(-5);
    N = 100;
    count = 0;
    % 定义函数原型，方便用于不同的方程，提高程序应用的普遍性
    syms f x;
    f = (x-1)^2*(2*x-1);
    % 求解方程  
    [answer,count] = advance_newton(x_initial,count,delta,f,N);
    fprintf("改进的牛顿法求非线性方程，初始值为 %f 解为 %f 迭代了 %d 次\n",x_initial,answer,count);
end
%%
% 方法五：拟牛顿法-秩1的拟牛顿法-逆Broyden法
% 题目：用拟牛顿法-逆Broyden求解下列非线性方程组的根，题目见实验报告册
% 注意，以下Fcn仅适用于3x3阶非线性方程组求解
% TODO 改成 NxN阶非线性方程组求解

if algorithm_index == 4
    % 定义迭代初始解向量
    X_initial = [1.0 1.0 1.0]'; 
    % 定义允许的误差以及最大的迭代次数
    delta = 0.5*10^(-5);
    N = 100;
    count = 0;
    % 定义函数原型，方便用于不同的方程，提高程序应用的普遍性
    syms f_1 f_2 f_3 x y z;
    f_1 = x*y-z^2-1;
    f_2 = x*y*z+y^2-x^2-2;
    f_3 = exp(x)+z-exp(y)-3;
    F = [f_1 f_2 f_3]';
    % 求系数矩阵A0
    A_initial = [diff(f_1,x) diff(f_1,y) diff(f_1,z); diff(f_2,x) diff(f_2,y) diff(f_2,z); diff(f_3,x) diff(f_3,y) diff(f_3,z)];
    x = X_initial(1);
    y = X_initial(2);
    z = X_initial(3);
    % 求系数矩阵H0
    H_initial = inv(eval(A_initial));
    % 求解方程  
    [answer,count] = quasi_newton(X_initial,H_initial,count,delta,F,N);
    fprintf("拟牛顿法求非线性方程组，初始值为 [%f %f %f]' 解为 [%f %f %f]' 迭代了 %d 次\n",X_initial(1),X_initial(2),X_initial(3),answer(1),answer(2),answer(3),count);
end
%% 
% 定义迭代函数实现二分法
function [answer,count] = dichotomy(next_x,next_y,count,delta,f)
answer = (next_x+next_y)/2;

x = next_y;
if eval(f) ==0
    answer = next_y;
elseif next_y - next_x > 2*delta
    count = count+1;
    x = answer;
    f1 = eval(f);
    x = next_x;
    f2 = eval(f);
    if f1*f2 > 0
        next_x = answer;
    else % 包含了端点值为解的情况
        next_y = answer;
    end
    [answer,count] = dichotomy(next_x,next_y,count,delta,f);
end
end
%% 
% 定义迭代函数实现牛顿法
function [answer,count] = newton(x_initial,count,delta,f,N)
answer = x_initial;
x = answer;
answer = answer - eval(f)/eval(diff(f));
count = count +1;
%if (abs(eval(f)) > delta)||(abs(eval(f)/eval(diff(f))) > delta)
if abs(eval(f)/eval(diff(f))) > delta
    if count < N
        [answer,count] = newton(answer,count,delta,f,N);
    else
        fprinf("Error, can not solve this equation in a limited count of %d",N);
    end
end
end
%%
% 定义迭代函数实现割线法
function [answer,count] = secant(x0_initial,x1_initial,count,delta,f,N)
answer = x1_initial;
answer_k = x1_initial;
answer_k_1 = x0_initial;
x = answer_k_1;
f0 = eval(f);
x = answer_k;
f1 = eval(f);
answer = answer_k - f1/(f1-f0)*(answer_k-answer_k_1);
count = count +1;
if abs(f1/(f1-f0)*(answer_k-answer_k_1)) > delta
    if count < N
        [answer,count] = secant(answer_k,answer,count,delta,f,N);
    else
        fprinf("Error, can not solve this equation in a limited count of %d",N);
    end
end
end
%% 
% 定义迭代函数实现改进的牛顿法
function [answer,count] = advance_newton(x_initial,count,delta,f,N)
answer = x_initial;
x = answer;
answer = answer - 2*eval(f)/eval(diff(f));
count = count +1;
fprintf("第%d次迭代，迭代解为%f 两次解的差值为%f delta为%f\n",count,answer,abs(2*eval(f)/eval(diff(f))),delta);
%if (abs(eval(f)) > delta)||(abs(eval(f)/eval(diff(f))) > delta)
if abs(2*eval(f)/eval(diff(f))) > delta
    if count < N
        [answer,count] = advance_newton(answer,count,delta,f,N);
    else
        fprintf("Error, can not solve this equation in a limited count of %d",N);
    end
end
end
%% 
% 定义迭代函数实现拟牛顿法-逆Broyden法
function [answer,count] = quasi_newton(X_initial,H_initial,count,delta,F,N)
answer = X_initial;
x = X_initial(1);
y = X_initial(2);
z = X_initial(3);
F_i = eval(F);
answer = answer - H_initial*F_i;
count = count + 1;
% 比较差值向量的 X(i+1) - X(i)的无穷范数与delta
if norm(H_initial*F_i,inf) > delta
    R =  - H_initial*eval(F);
    x = answer(1);
    y = answer(2);
    z = answer(3);
    F_i_1 = eval(F);
    Y = F_i_1 - F_i;
    H_initial = H_initial + (R-H_initial*Y)*(R'*H_initial)/(R'*H_initial*Y);
    if count < N
        [answer,count] = quasi_newton(answer,H_initial,count,delta,F,N);
    else
        fprinf("Error, can not solve this equation set in a limited count of %d",N);
    end
end
end

%% ------------------END OF THE FILE------------------