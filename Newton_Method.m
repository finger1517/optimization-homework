function [xmin,fmin]=Newton_Method(F,gradF,HessianF,x0,TolGrad,MaxIter)
%
%  Input
%  F........ function handle for objective function F(x) with input argument, x
%  gradF... function handle for gradient of objective function with input argument, x
%  x0....... design variable vector initial guess for starting point
%  TolGrad.. Tolerance for norm of objective gradient to be zero
%  MaxIter.. Maximum number of iterations
%
%  Output
%  x........ final design variable vector found to minimize objective function, F(x)
%  f1....... final objective function minimum value

if nargin<5 || isempty(TolGrad), TolGrad=5e-4; end
if nargin<6 || isempty(MaxIter), MaxIter=100;  end

% using Newton method

%% Initialize loop parameters
iter = 0;
f0 = F(x0);
c0 = gradF(x0);
H = HessianF(x0);
c  = c0;
Converged = norm(c) < TolGrad;
alpha = 1;

disp('iter ----   alpha ----   f(x)  ----   norm(gradient)');
fprintf('%4.0f  ----  %6.6f ----   %8.4f ----   %8.4f\n',[iter, 0, f0, norm(c)]);

%% Search direction and line search iterations

while iter<MaxIter && ~Converged

    iter = iter + 1;
%     fprintf('-----------------current iteration is %d---------------------\n',iter)
    
	d = -inv(H)*c;

        alpha = linesearch(F,gradF,x0,d);
        f1 = F(x0);
	x = x0 + alpha*d;
    
   %renew the current point 
   x = x0 + alpha*d;
   c = gradF(x);
   H = HessianF(x);
   Converged = norm(c) < TolGrad;
   x0 = x;
   c0 = c;
   
   fprintf('%4.0f---- %6.4f ----%8.4f---- %8.4f\n',[iter, alpha, f1, norm(c)])
%    if ~Converged
%        disp('not converged, current norm of gradient equals to:');
%        disp(norm(c));
%    end
   
%    disp('is iter<Maxiter?');
%    disp(iter<MaxIter);
   
end


xmin = x0;
fmin = F(x0);
disp('---------------------------');
disp('the solution:');
disp(x0);
disp('minimal value:');
disp(fmin);
end
