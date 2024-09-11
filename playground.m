clear all
x = linspace(0,2,6);
y = sin(pi*x);

k = 6; 
sp = spapi( optknt(x,k), x, y );

xx = linspace(0,2,1001); % desired sampling
y_estimated = fnval(xx,sp)';

y_estimated = y_estimated';
y_actual = sin(pi * xx);
plot (x,y,'*',xx,y_estimated, xx, y_actual)

sqrth = 1.0/sqrt(6); % error factor for spatial l2 norm
errinf = norm((y_estimated-y_actual),inf);     % estimate of infinity error 
err2   = sqrth*norm((y_estimated-y_actual),2); % estimate of 2 norm

fprintf('%s | n = %3i |  inf error = %8.2e | 2 norm error = %8.2e \n','AA', 6,errinf,err2);

y_estimated(1:5)
y_actual(1:5)