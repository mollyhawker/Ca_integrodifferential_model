function [Y_new] = RK4(dydt, t, Y, dt)
% RK4 Implementation of 4th-order Runge Kutta
   
    YK1=dt*dydt(t, Y);
    
    Y1=Y+YK1/2;
    
    YK2=dt*dydt(t, Y1);
    
    Y2=Y+YK2/2;
    
    YK3=dt*dydt(t, Y2);
    
    Y3=Y+YK3;
    
    YK4=dt*dydt(t,Y3);
    
    Y_new=Y+(YK1+2*YK2+2*YK3+YK4)/6;
    
end

