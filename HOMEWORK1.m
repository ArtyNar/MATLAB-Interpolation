%clear all
%close all      
warning('off') 

function [] = do_plotting(yplot, x, expx)
   % Creates eveny spaced points for error sampling
   xplot = linspace(0,2,1001);  
   % Calculates the true values for the function for error sampling
   expplot = exp(xplot);
   
   % Setting up a figure with two graphs
   figure

   % Plots x mesh, true values of y
   subplot(2,1,1);
   plot(xplot, expplot, 'r');

   % Adds legend to the plot
   title('True graph')
   xlabel(' x ');
   ylabel('exp(x)');

   % Plots x mesh, approximated values of y, and approximation points
   subplot(2,1,2);
   plot(xplot, yplot, x, expx, 'o');
   
   % Adds legend to the plot
   title('Interpolated graph')
   legend('Interpolated graph', 'Interpolation points');
   xlabel(' x ');
   ylabel('exp(x) - Approximated');
end

% A function for Polynomial Interpolation using Vandermonde matrix 
% with n evenly spacedpoints
function [] = evenly_spaced_vander(n)
   % Creates eveny spaced points for error sampling
   xplot = linspace(0,2,1001);  
   % Calculates the true values for the function for error sampling
   expplot = exp(xplot);
   % Creates a mesh of x for the n points
   x = linspace(0,2,n)';

   % Caclulates data values at x 
   expx = exp(x);         

   % Creates a Vandermonde marix V
   V = fliplr(vander(x));           

   % Solves for coefficients of powers  of x
   a  = V\expx;                     

   % Forming the polynomial
   yplot = ones(1,1001)*a(n); 
   for j = 1:n-1
      yplot = yplot.*xplot + a(n-j); % Note vector operations of Horners scheme  
   end

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('Vander (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);

   % Uses a function to simplify plotting
   do_plotting(yplot, x, expx);
end

function [] = chebyshev_spaced_vander(n)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);
   

   % Defines the Chebyshev data points on range [0-2]
   x = zeros(n,1);
   for k = 1:n
      x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
   end

   % Calculates the values at the Chebyshev points
   expx = exp(x);

   % Creates a Vandermonde marix V
   V = fliplr(vander(x));           

   % Solves for coefficients of powers  of x
   a  = V\expx;                     

   % Forming polynomial
   yplot = ones(1,1001)*a(n); 
   for j = 1:n-1
      yplot = yplot.*xplot + a(n-j);  % note vector operations of Horners scheme
   end

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('Vander (Chebyshev) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);

   % Uses a function to simplify plotting
   do_plotting(xplot, expplot, yplot, x, expx);
end

function [] = evenly_spaced_bary(n)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Creates a mesh of x for the n evenly-spaced points
   x = linspace(0,2,n)';

   % Caclulates data values at x 
   expx = exp(x);     
   
   vk = baryWeights(x); % Barycentric weights
   yplot = bary(xplot, expx, x, vk);

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('Barycentric (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);

   do_plotting(yplot, x, expx);
end

function [] = chebyshev_spaced_bary(n)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Defines the Chebyshev data points on range [0-2]
   x = zeros(n,1);
   for k = 1:n
      x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
   end

   % Caclulates data values at x 
   expx = exp(x);     
   
   vk = baryWeights(x); % Barycentric weights
   yplot = bary(xplot, expx, x, vk);

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('Barycentric (chebyshev-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);

   do_plotting(yplot, x, expx);
end

function [] = evenly_spline(n, k)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Creates a mesh of x for the n evenly-spaced points
   x = linspace(0,2,n);

   % Caclulates data values at x 
   expx = exp(x);

   sp = spapi(optknt(x,k), x, expx);

   yplot = fnval(xplot,sp)';

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('%d-Order Spline (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n',k, n,errinf,err2);
   
   do_plotting(yplot, x, expx);
end

function [] = chebyshev_spline(n, k)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Defines the Chebyshev data points on range [0-2]
   x = zeros(n,1);
   for k = 1:n
      x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
   end

   x = x';
   % Caclulates data values at x 
   expx = exp(x);

   sp = spapi(optknt(x,k), x, expx);

   yplot = fnval(xplot,sp)';

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('%d-Order Spline (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n',k, n,errinf,err2);
   
   do_plotting(yplot, x, expx);
end

function [] = evenly_pchip(n)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Creates a mesh of x for the n evenly-spaced points
   x = linspace(0,2,n)';

   % Caclulates data values at x 
   expx = exp(x);

   % Values using cubic spline
   yplot = interp1(x,expx, xplot,'pchip');

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('Pchip (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);
   
   do_plotting(yplot, x, expx);
end

function [] = chebyshev_pchip(n)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Defines the Chebyshev data points on range [0-2]
   x = zeros(n,1);
   for k = 1:n
      x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
   end

   % Caclulates data values at x 
   expx = exp(x);

   % Values using cubic spline
   yplot = interp1(x,expx, xplot,'pchip');

   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
   fprintf('Pchip (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);
   
   do_plotting(yplot, x, expx);
end

fprintf('=========================================\n')

% for n = [6,11,21,41,81,161,321,641]
%    fprintf('Vandermonde matrix with %d evenly-spaced points:\n', n);
%    evenly_spaced_vander(n)
% end

%  fprintf('=========================================\n')

% for n = [6,11,21,41,81,161,321,641]
%    fprintf('Vandermonde matrix with %d Chebyshev points:\n', n);
%    chebyshev_spaced_vander(n)
% end


% A function for Polynomial Interpolation using Vandermonde matrix 
% with n evenly spacedpoints
function [] = trystuff()
   % Creates eveny spaced points for error sampling
   xplot = linspace(0,2,1001);  
   % Calculates the true values for the function for error sampling
   expplot = exp(xplot);
   
   % Setting up a figure for future graphing
   f = figure;
   f.Name = 'Interpolation';
   f.Position(1:4) = [200 200 900 600];
   % Plots x mesh, true values of y
   subplot(3,3,1);
   plot(xplot, expplot, 'r');
   i = 1;

   for n = [6,11,21,41,81,161,321,641]
      % Creates a mesh of x for the n points
      x = linspace(0,2,n)';

      % Caclulates data values at x 
      expx = exp(x);         

      % Creates a Vandermonde marix V
      V = fliplr(vander(x));           

      % Solves for coefficients of powers  of x
      a  = V\expx;                     

      % Forming the polynomial
      yplot = ones(1,1001)*a(n); 
      for j = 1:n-1
         yplot = yplot.*xplot + a(n-j); % Note vector operations of Horners scheme  
      end

      % Estimating errors
      sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
      errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
      err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm
      fprintf('Vander (evenly-spaced) n = %i  inf error= %8.2e 2 norm error = %8.2e \n', n,errinf,err2);


      % Adds legend to the plot
      title('True graph')
      xlabel(' x ');
      ylabel('exp(x)');

      % Plots x mesh, approximated values of y, and approximation points
      subplot(3,3,i+1);
      plot(xplot, yplot, x, expx, '*');
      
      % Adds legend to the plot
      title('Interpolated graph')
      legend('Interpolated graph', 'Interpolation points');
      xlabel(' x ');
      ylabel('exp(x) - Approximated');

      i = i + 1;
   end
end

trystuff()