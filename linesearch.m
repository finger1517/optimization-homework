function alpha = linesearch(F,gradF,x,d)
%         using amijio condition to do the linesearch
        alpha = 1;
        c = 0.01;
        tao = 0.5;
        while 1
                fn = F(x + alpha*d);
                if fn > F(x)-alpha*c*d'*d
                       alpha = tao * alpha;
                else
%                         alpha = max(0.05, alpha);
                        break;
                end
        end
                
                
end
        