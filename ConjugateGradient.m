function [xmin,fmin]=ConjugateGradient(F,gradF,x0,TolGrad,MaxIter)
%
%  Input
%  F........ function handle for objective function F(x) with input argument, x
%  gradF... function handle for gradient of objective function with input argument, x
%  x0....... design variable vector initial guess for starting point
%  TolGrad.. Tolerance for norm of objective gradient to be zero
%  MaxIter.. Maximum number of iterations
%
%  Output
%  xmin........ final design variable vector found to minimize objective function, F(x)
%  fmin....... final objective function minimum value

if nargin<4 || isempty(TolGrad), TolGrad=1e-3; end
if nargin<5 || isempty(MaxIter), MaxIter=20; end

%% Change following steepest descent algorithm to Conjugate Gradient Method
% initialize loop parameters
iter = 0;
x = x0;
f0   = F(x0);
g0    = gradF(x0);
g = g0;
g_old = g;
Converged = norm(g0) < TolGrad;
disp('   iter --------------alpha-------------- f(x)--------------  norm(g)--------------')
fprintf(' %6.0f-------------- %5.3f-------------- %8.4f-------------- %8.4f\n',[iter, 0, f0, norm(g0)])
alpha = 1;


while iter<MaxIter && ~Converged
        if iter == 0
                p = -g0;
        else
                beta = g'*(g-g_old)/(g_old'*g_old);
                p = -g + beta*p;
        end
        alpha = linesearch(F, gradF, x, p);
        x = x + alpha * p;
        g_old = g;
        g = gradF(x);
        Converged = norm(g)<TolGrad;
        
        iter = iter + 1;
        f = F(x);
        disp('   iter --------------alpha----------------- f(x)-------------------norm(g)')
        fprintf(' %6.0f-------------- %5.3f-------------- %8.4f-------------- %8.4f\n',[iter, alpha, f, norm(g)])
end
xmin = x;
disp('xmin=');
disp(xmin);
fmin = F(xmin);
end