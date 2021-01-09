
%% 
% 哈工大数值分析2020年秋研究生，上机实验
% 第四部分 | 龙贝格积分法/Romberg Integral
% 时间: 2020/10/29
% 学生: 20S****** ***
% ----------------------------------------------------------
% 1、【龙贝格积分法计算以下积分的近似值】

%%
% 定义被积函数原型及精度要求
syms f x;
delta = 7*10^(-6);
T = zeros(8);
% 求解定积分1
f = x^3;
a = 6; b = 100;m = 1;
% f = 4/(1+x^2);
% a = 0; b = 1;m = 1;
[answer,m] = Romberg(f,a,b,m,delta,T);
fprintf("\n求解定积分1，解为 %6f    m =  %d  验证解为%f\n-------------------------------------\n",answer, m, int(eval(f),x,a,b));
% 求解定积分2
f = sin(x)/x;
a = 0; b = 1;m = 1;
[answer,m] = Romberg(f,a,b,m,delta,T);
fprintf("\n求解定积分2，解为 %f    m =  %d  验证解为%f\n-------------------------------------\n",answer, m, int(eval(f),x,a,b));
% 求解定积分3
f = sin(x^2);
a = 0; b = 1;m = 1;
[answer,m] = Romberg(f,a,b,m,delta,T);
fprintf("\n求解定积分3，解为 %f    m =  %d  验证解为%f\n-------------------------------------\n",answer, m, int(eval(f),x,a,b));
%%
% 定义Romberg积分迭代函数
function [answer,m] = Romberg(f,a,b,m,delta,T)
if m == 1
    x = a;
    if x == 0
        syms x;
        fa = limit(eval(f),x,0);
    else
        fa = eval(f); 
    end
    x = b;
    fb = eval(f);
    T(1,1) = 0.5*(b-a)*(fa+fb);
end
for i = 0:m
    % 计算龙贝格T-数表
    if  i == 0 % T型求积公式
        sum = 0;
        for j = 0:2^(m-1)-1
            x = a + (b-a)/2^(m-1)*(j+0.5);
            if x == 0
                syms x;
                sum = sum + limit(eval(f),x,0);
            else
                sum = sum + eval(f);
            end
        end
        T(m+1,i+1) = 0.5*T(m,i+1)+0.5*(b-a)/2^(m-1)*sum;
    else % 高阶求积公式
        T(m+1,i+1) = (4^i*T(m+1,i)-T(m,i))/(4^i-1);
    end
end
% 计算龙贝格数表上的对角线最后两元素之差，与给定的精度进行比较
if abs(T(m+1,m+1)-T(m,m)) <= delta
    answer = T(m+1,m+1);
    fprintf("求解完成!  T-数表为\n\n");
    for i = 0:m
        for j = 0:i
            fprintf("%f  ",T(i+1,j+1));
        end
        fprintf("\n");
    end
else
    % 未达到精度要求，继续迭代求解
    [answer,m] = Romberg(f,a,b,m+1,delta,T);
end
end
%% ------------------END OF THE FILE------------------