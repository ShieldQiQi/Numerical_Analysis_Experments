
%% 
% ��������ֵ����2020�����о������ϻ�ʵ��
% ��һ���� | �����Է��������
% ʱ��: 2020/10/15
% ѧ��: 20S****** ***
% ----------------------------------------------------------
% 1�������ַ���
% 2����ţ�ٷ���
% 3�������߷���
% 4�����Ľ���ţ�ٷ���
% 5������ţ�ٷ���
%%
% define the algorithm for which be used to solve the nolinear equation
% where the variable can be 0��1��2��3��4 , each one corresponds to the algorithm above

algorithm_index = 4;

%%
% ����һ�����ַ�
% ��Ŀ���ö��ַ����㷽�� sin(x)-pow(x,2)/2 = 0 ��(1,2)�ڵĸ��Ľ���ֵ��Ҫ���=0.5*10^(-5)
    
if algorithm_index == 0
    % ��֪����equation(x)������(1,2)�������ɵ����ҷ��̵ĸ�������(1,2)�ڴ�����Ψһ
    a = 1;
    b = 2;
    count = 0;
    % ������������
    delta = 0.5*10^(-5);
    % ���庯��ԭ�ͣ��������ڲ�ͬ�ķ��̣���߳���Ӧ�õ��ձ���
    syms f x;
    f = sin(x)-x^2/2;
    % ���������Լ���������
    [answer,count] = dichotomy(a,b,count,delta,f);
    fprintf("���ַ�������Է��̣���Ϊ %f ������ %d ��\n",answer,count);
end
%%
% ��������ţ�ٷ�
% ��Ŀ����ţ�ٷ�������з����Է��̵ĸ�����Ŀ��ʵ�鱨���

if algorithm_index == 1
    % ����3�����̵ĳ�ʼֵ
    x1_initial = 0.5;
    x2_initial = 1;
    x3_initial_1 = 0.45;
    x3_initial_2 = 0.65;
    % �������������Լ����ĵ�������
    delta = 0.5*10^(-5);
    N = 100;
    % ���庯��ԭ�ͣ��������ڲ�ͬ�ķ��̣���߳���Ӧ�õ��ձ���
    syms f x;
    % ��ⷽ��1
    count = 0;
    f = x*exp(x)-1;
    [answer,count] = newton(x1_initial,count,delta,f,N);
    fprintf("ţ�ٷ�������Է���1����ʼֵΪ %f ��Ϊ %f ������ %d ��\n",x1_initial,answer,count);
    % ��ⷽ��2
    count = 0;
    f = x^3-x-1;
    [answer,count] = newton(x2_initial,count,delta,f,N);
    fprintf("ţ�ٷ�������Է���2����ʼֵΪ %f ��Ϊ %f ������ %d ��\n",x2_initial,answer,count);
    % ��ⷽ��3
    count = 0;
    f = (x-1)^2*(2*x-1);
    [answer,count] = newton(x3_initial_1,count,delta,f,N);
    fprintf("ţ�ٷ�������Է���3����ʼֵΪ %f ��Ϊ %f ������ %d ��\n",x3_initial_1,answer,count);
    count = 0;
    [answer,count] = newton(x3_initial_2,count,delta,f,N);
    fprintf("ţ�ٷ�������Է���3����ʼֵΪ %f ��Ϊ %f ������ %d ��\n",x3_initial_2,answer,count);
end
%%
% �����������߷�/��������
% ��Ŀ���ø��߷�������з����Է��̵ĸ�����Ŀ��ʵ�鱨���

if algorithm_index == 2
    % ���������ʼ��
    x0_initial = 0.4; 
    x1_initial = 0.6; 
    % �������������Լ����ĵ�������
    delta = 0.5*10^(-5);
    N = 100;
    count = 0;
    % ���庯��ԭ�ͣ��������ڲ�ͬ�ķ��̣���߳���Ӧ�õ��ձ���
    syms f x;
    f = x*exp(x)-1;
    % ��ⷽ��  
    [answer,count] = secant(x0_initial,x1_initial,count,delta,f,N);
    fprintf("���߷�������Է��̣���ʼֵx0Ϊ %f ��ʼֵx1Ϊ %f ��Ϊ %f ������ %d ��\n",x0_initial,x1_initial,answer,count);
end
%%
% �����ģ��Ľ���ţ�ٷ�
% ��Ŀ���øĽ���ţ�ٷ�������з����Է��̵ĸ�����Ŀ��ʵ�鱨���

if algorithm_index == 3
    % ���������ʼ��
    x_initial = 0.55; 
    % �������������Լ����ĵ�������
    delta = 0.5*10^(-5);
    N = 100;
    count = 0;
    % ���庯��ԭ�ͣ��������ڲ�ͬ�ķ��̣���߳���Ӧ�õ��ձ���
    syms f x;
    f = (x-1)^2*(2*x-1);
    % ��ⷽ��  
    [answer,count] = advance_newton(x_initial,count,delta,f,N);
    fprintf("�Ľ���ţ�ٷ�������Է��̣���ʼֵΪ %f ��Ϊ %f ������ %d ��\n",x_initial,answer,count);
end
%%
% �����壺��ţ�ٷ�-��1����ţ�ٷ�-��Broyden��
% ��Ŀ������ţ�ٷ�-��Broyden������з����Է�����ĸ�����Ŀ��ʵ�鱨���
% ע�⣬����Fcn��������3x3�׷����Է��������
% TODO �ĳ� NxN�׷����Է��������

if algorithm_index == 4
    % ���������ʼ������
    X_initial = [1.0 1.0 1.0]'; 
    % �������������Լ����ĵ�������
    delta = 0.5*10^(-5);
    N = 100;
    count = 0;
    % ���庯��ԭ�ͣ��������ڲ�ͬ�ķ��̣���߳���Ӧ�õ��ձ���
    syms f_1 f_2 f_3 x y z;
    f_1 = x*y-z^2-1;
    f_2 = x*y*z+y^2-x^2-2;
    f_3 = exp(x)+z-exp(y)-3;
    F = [f_1 f_2 f_3]';
    % ��ϵ������A0
    A_initial = [diff(f_1,x) diff(f_1,y) diff(f_1,z); diff(f_2,x) diff(f_2,y) diff(f_2,z); diff(f_3,x) diff(f_3,y) diff(f_3,z)];
    x = X_initial(1);
    y = X_initial(2);
    z = X_initial(3);
    % ��ϵ������H0
    H_initial = inv(eval(A_initial));
    % ��ⷽ��  
    [answer,count] = quasi_newton(X_initial,H_initial,count,delta,F,N);
    fprintf("��ţ�ٷ�������Է����飬��ʼֵΪ [%f %f %f]' ��Ϊ [%f %f %f]' ������ %d ��\n",X_initial(1),X_initial(2),X_initial(3),answer(1),answer(2),answer(3),count);
end
%% 
% �����������ʵ�ֶ��ַ�
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
    else % �����˶˵�ֵΪ������
        next_y = answer;
    end
    [answer,count] = dichotomy(next_x,next_y,count,delta,f);
end
end
%% 
% �����������ʵ��ţ�ٷ�
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
% �����������ʵ�ָ��߷�
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
% �����������ʵ�ָĽ���ţ�ٷ�
function [answer,count] = advance_newton(x_initial,count,delta,f,N)
answer = x_initial;
x = answer;
answer = answer - 2*eval(f)/eval(diff(f));
count = count +1;
fprintf("��%d�ε�����������Ϊ%f ���ν�Ĳ�ֵΪ%f deltaΪ%f\n",count,answer,abs(2*eval(f)/eval(diff(f))),delta);
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
% �����������ʵ����ţ�ٷ�-��Broyden��
function [answer,count] = quasi_newton(X_initial,H_initial,count,delta,F,N)
answer = X_initial;
x = X_initial(1);
y = X_initial(2);
z = X_initial(3);
F_i = eval(F);
answer = answer - H_initial*F_i;
count = count + 1;
% �Ƚϲ�ֵ������ X(i+1) - X(i)���������delta
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