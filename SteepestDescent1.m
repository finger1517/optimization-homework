function [xmin,fmin]=SteepestDescent1(F,gradF,x0,TolGrad,MaxIter)
%
%  Input
%  F........ function handle for objective function F(x) with input argument, x
%  gradF... function handle for gradient of objective function with input argument, x
%  x0....... design variable vector initial guess for starting point
%  TolGrad.. Tolerance for norm of objective gradient to be zero
%  MaxIter.. Maximum number of iterations
%  Output
%  xmin........ final design variable vector found to minimize objective function, F(x)
%  fmin....... final objective function minimum value

if nargin<4 || isempty(TolGrad), TolGrad=1e-2; end
if nargin<5 || isempty(MaxIter), MaxIter=20; end

iter = 0;
p = gradF(x0);
x = x0;
converged = norm(p) < TolGrad;
% alpha = 1e-6*norm(p);

while iter < MaxIter &&  ~converged 
        iter = iter +1;
        d = -p;
%        alpha = linesearch(F, gradF, x, d);
         alpha = simplelinesearch(F, gradF, x, p);
% too slow         alpha = AlphaHelper(F, x, d);
%         f = @(alpha) F(x-alpha*p);
%         alphaUpper = bracket( f, 0, 0.1*alpha );
%         [alpha,~] = fminbnd( f, 0, alphaUpper );
        x = x - alpha * p;
        
        f = F(x);
        p = gradF(x);
        converged = norm(p) < TolGrad;
        fprintf('--------------iter:%d-----------------norm(g):%f-------------------value:%f--------alpha:%f\n',iter,norm(p),f,alpha)
end

disp(x);
xmin = x;
fmin = F(x);
end



