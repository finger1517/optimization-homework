function [x,fmin]=BFGS_Quasi_Newton(F,gradF,x0,TolGrad,MaxIter)
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

if nargin<4 || isempty(TolGrad), TolGrad=1e-3; end
if nargin<5 || isempty(MaxIter), MaxIter=100;  end

% Change following steepest descent algorithm to Quasi-Newton Method 
% using BFGS Hessian updates.

%% Initialize loop parameters
iter = 0;
dim_x = length(x0);
x = x0;
grad1 = gradF(x0);
B = eye(dim_x);
Converged = norm(grad1) < TolGrad;
f0 = F(x0);
alpha = 1;
Fconverge = false;
Ftol = 1e-7;

%% Search direction and line search iterations
fprintf('iter: %d  --------------   f(x): %f   ------------------|gradient|:%f   ---------------alpha:%f\n',iter,f0,norm(grad1),alpha);
while iter<MaxIter && ~Converged && ~Fconverge
	iter = iter + 1;
%         
%         disp('x=');
%         disp(x);
%         disp(norm(grad1));
%         disp(F(x));
        s = -1 * B *grad1;
        s = s/norm(s);
        
        alpha = linesearch(F,gradF,x,s);
        d = alpha*s;
        fold = F(x);
        x  = x + d;
        f = F(x);
        grad2 = gradF(x);
        y = grad2 - grad1;
        grad1 = grad2;
        
        %using BFGS Method to update the H
        if (y'*s > 0)
                B = BFGS_update(B, s, y);
        end
        
        delta = abs(f - fold);
        
        Converged = norm(grad1) < TolGrad;
        Fconverge = delta < Ftol;
        
%         if Converged
%                 disp('converged');
%         end
%         
        fprintf('iter:%d   --------------f(x):   %f   ------------------|gradient|:%f   ---------------alpha:%f\n',iter,F(x),norm(grad1),alpha);
end

xmin = x;
fmin = F(x);
% disp('xmin=');
% disp(xmin);
% fprintf('%4.0f %6.4f %8.4f %8.4f\n',[iter, alpha, fmin, norm(c)])
end


function [B] = BFGS_update(H, s, y)
% this function is to update the inverted function B_k+1
% H  = inv(B_k) is the last inverted Hessian matrix
        N = length(s);
        a = s*y';
        b = y'*s;
        c = y*s';
        d = s*s';
        B = (eye(N) - a/b)' * H * (eye(N) - c/b) + d/b;
end

function alpha = linesearch(F,gradF,x,d)
%         using amijio condition to do the linesearch
        alpha = 1;
        iter = 0;
        c = 0.001;
        tao = 0.5;
        while 1
                fn = F(x + alpha*d);
                if fn > F(x)+alpha*c*d'*gradF(x)
                       alpha = tao * alpha;
                else
%                         alpha = max(0.0001, alpha);
                        break;
                end
        end
                
                
end
        
        

    