function CO=Test_opt4(x)
global t y CO

CO=0;Ori=0;
if length(x)==3
    f=x(1)*sin(x(2)*pi*t)+x(3);
    for i=1:length(t)        
        CO=CO+(y(i)-f(i))^2;  
    end
    
elseif length(x)==4
    f=x(1)*sin(x(2)*pi*t)+x(3)+x(4)*t;
    for i=1:length(t)
        CO=CO+abs(y(i)-f(i));       
    end
   
elseif length(x)==5
    f=x(1)*sin(x(2)*pi*t-x(5))+x(3)+x(4)*t;
    for i=1:length(t)
        CO=CO+abs(y(i)-f(i));
        Ori=Ori+abs(y(i));        
    end
    CO=CO/Ori;
end

% xx=x;
end