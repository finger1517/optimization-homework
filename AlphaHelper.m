function alpha_opt = AlphaHelper(F, x, d)
    a = 0;
    b = 10;
    eps = 0.01;
    
    %initial
    f = F(x);
%     d = -H * g;
%     d = d/norm(d);
    
    while ((b-a)>eps)
        xL = a + (b-a)/4;
        xR = b - (b-a)/4;
        
%         fprintf('a and b = %f, %f \n',a,b);
        
        L = x + xL * d;
        R = x + xR * d;
        
%         fprintf('L and R = %f, %f \n',L,R);
        
%         disp('FL and FR');
%         disp(F(L));
%         disp(F(R));
        if(F(L)>F(R))
            a = xL;
        else
            b = xR;
        end
    end
    
    alpha_opt = (a+b)/2;
    
end