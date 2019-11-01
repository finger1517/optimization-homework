%% 初始化
clc
clear
% Rosenbrock's banana function
x0 = zeros(10,1);
x1 = ones(10,1);
x2 = 1.5*x1;
F = @(x)rosen(x);
gradF = @(x)grad(x);
hessianF = @(x)HessianMatrix(x);



%% 测试
[xmin,fmin]=ConjugateGradient(F,gradF,x0,0.01,1000);%525 steps
% [xmin,fmin]=SteepestDescent1(F,gradF,x0, 0.01,5000); %  4241 steps,和linesearch的参数tao有关
%  [xmin,fmin]=Newton_Method(F,gradF,hessianF,x0,0.01,25); % 24 steps
%[xmin,fmin]=BFGS_Quasi_Newton(F,gradF,x0,0.01,10000) %1098 steps
