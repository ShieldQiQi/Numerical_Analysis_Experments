
%% 
% ��������ֵ����2020�����о������ϻ�ʵ��
% ���Ĳ��� | ��������ַ�/Romberg Integral
% ʱ��: 2020/10/29
% ѧ��: 20S****** ***
% ----------------------------------------------------------
% 1������������ַ��������»��ֵĽ���ֵ��

%%
% ���屻������ԭ�ͼ�����Ҫ��
syms f x;
delta = 7*10^(-6);
T = zeros(8);
% ��ⶨ����1
f = x^3;
a = 6; b = 100;m = 1;
% f = 4/(1+x^2);
% a = 0; b = 1;m = 1;
[answer,m] = Romberg(f,a,b,m,delta,T);
fprintf("\n��ⶨ����1����Ϊ %6f    m =  %d  ��֤��Ϊ%f\n-------------------------------------\n",answer, m, int(eval(f),x,a,b));
% ��ⶨ����2
f = sin(x)/x;
a = 0; b = 1;m = 1;
[answer,m] = Romberg(f,a,b,m,delta,T);
fprintf("\n��ⶨ����2����Ϊ %f    m =  %d  ��֤��Ϊ%f\n-------------------------------------\n",answer, m, int(eval(f),x,a,b));
% ��ⶨ����3
f = sin(x^2);
a = 0; b = 1;m = 1;
[answer,m] = Romberg(f,a,b,m,delta,T);
fprintf("\n��ⶨ����3����Ϊ %f    m =  %d  ��֤��Ϊ%f\n-------------------------------------\n",answer, m, int(eval(f),x,a,b));
%%
% ����Romberg���ֵ�������
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
    % ����������T-����
    if  i == 0 % T�������ʽ
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
    else % �߽������ʽ
        T(m+1,i+1) = (4^i*T(m+1,i)-T(m,i))/(4^i-1);
    end
end
% ���������������ϵĶԽ��������Ԫ��֮�������ľ��Ƚ��бȽ�
if abs(T(m+1,m+1)-T(m,m)) <= delta
    answer = T(m+1,m+1);
    fprintf("������!  T-����Ϊ\n\n");
    for i = 0:m
        for j = 0:i
            fprintf("%f  ",T(i+1,j+1));
        end
        fprintf("\n");
    end
else
    % δ�ﵽ����Ҫ�󣬼����������
    [answer,m] = Romberg(f,a,b,m+1,delta,T);
end
end
%% ------------------END OF THE FILE------------------