function alpha = simplelinesearch(F,gradF,x,p)
        beta = 0.1;
        alpha = 10;
        while 1
                if F(x - alpha * p) > F(x)
                        alpha = beta * alpha;
                else
%                         alpha = max(alpha, 0.01);
                        break;
                end
        end
end