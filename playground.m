for i = [6,11,21,41,81,161,321,641]
n = i ; % = 2^(i-2)*10+1;

x = linspace(0,1,n)';
sinx = sin(pi*x);
A = fliplr(vander(x));
a = A\sinx;
nplot = 10001; sqrth = 1.0/sqrt(n);
xplot = linspace(0,1,nplot);
sinplot= sin(pi*xplot);
yplot = ones(1,nplot)*a(n);
for j = 1:n-1
yplot = yplot.*xplot + a(n-j); % note vector operations
end
errinf = norm((yplot-sinplot),inf);
err2 = sqrth*norm((yplot-sinplot),2);
fprintf(' n= %i infinity error= %8.2e 2 norm error =%8.2e \n', n,errinf,err2);
end
