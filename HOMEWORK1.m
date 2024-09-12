%clear all
%close all      
warning('off') 

% A helper function that estimates errors
%   algorithm_name:  the text to display in printf
%   n:               number of points used
%   yplot:           estimated values
%   expplot:         true values
function [] = display_errors(algorithm_name, n,yplot, expplot)
   % Estimating errors
   sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
   errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
   err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm

   fprintf('%s | n = %3i |  inf error = %8.2e | 2 norm error = %8.2e \n',algorithm_name, n,errinf,err2);
end

% A function that can be used to display a single graph
%     yplot:   estimated values
%     x:       set of x values at points used for estimation 
%     y:       set of y values at points used for estimation
function [] = do_plotting_singular(yplot, x, expx)
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

% A function used to display 8 graphs at the same time
% It accepts a single set of values x and yplot, but when 
% used in a loop all the graphs will get fused in a figure
%     i:     used for positioning of the graph in a figure
%     n:     number of points used to interpolate
%     yplot: estimated y-valies
%     x:     points intialy used to do interpolation
function [] = do_plotting_all(i, n, yplot, x)
   % Creates eveny spaced points for error sampling
   xplot = linspace(0,2,1001);  

   % Caclulates data values at x 
   expx = exp(x);  

   % Plots x mesh, approximated values of y, and approximation points
   subplot(2,4,i);
   plot(xplot, yplot, x, expx, 'o', 'MarkerSize',3);
   
   % Adds legend to the plot
   title(strcat('Interpolated with n= ', num2str(n)))
   legend('Interpolated graph', 'Interpolation points');
   xlabel(' x ');
   ylabel('exp(x) - Approximated');
end

% A function for Polynomial Interpolation using Vandermonde matrix 
% with n evenly spaced points
function [yplot] = evenly_spaced_vander(n, x)
   % Creates eveny spaced points for error sampling
   xplot = linspace(0,2,1001);  
   % Calculates the true values for the function for error sampling
   expplot = exp(xplot);

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

   % The function will calculate and display errors
   display_errors('Vandermonde (evenly)',n, yplot, expplot);

   % Uses a function to simplify plotting
   % do_plotting_singular(yplot, x, expx);

end

% A function for Polynomial Interpolation using Vandermonde matrix 
% with n Chebyshev-spaced points
function [yplot] = chebyshev_spaced_vander(n, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

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

   % The function will calculate and display errors
   display_errors('Vandermonde (Chebyshev)',n, yplot, expplot);

   % Uses a function to simplify plotting
   % do_plotting_singular(yplot, x, expx);
end

% A function for Polynomial Lagrange Interpolation using bary/baryweights routines 
% with n evenly spaced points
function [yplot] = evenly_spaced_bary(n, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Caclulates data values at x 
   expx = exp(x);     
   
   vk = baryWeights(x); % Barycentric weights
   yplot = bary(xplot, expx, x, vk);

   % The function will calculate and display errors
   display_errors('Barycentric (evenly)',n, yplot, expplot);

   % do_plotting_singular(yplot, x, expx);
end

% A function for Polynomial Lagrange Interpolation using bary/baryweights routines 
% with n Chebyshev spaced points
function [yplot] = chebyshev_spaced_bary(n, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Caclulates data values at x 
   expx = exp(x);     
   
   vk = baryWeights(x); % Barycentric weights
   yplot = bary(xplot, expx, x, vk);

   % The function will calculate and display errors
   display_errors('Barycentric (Chebyshev)',n, yplot, expplot);

   % do_plotting_singular(yplot, x, expx);
end

% A function for spline interpolation 
% with n evenly spaced points
%      k-order pline
function [yplot] = evenly_spline(n, k, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   x = x';
   % Caclulates data values at x 
   expx = exp(x);

   sp = spapi(optknt(x,k), x, expx);

   yplot = fnval(xplot,sp);

   % The function will calculate and display errors
   display_errors('Spline (evenly)',n, yplot, expplot);
   
   %do_plotting_singular(yplot, x, expx);
end

% A function for spline interpolation 
% with n Chebyshev spaced points
%      k-order pline
function [yplot] = chebyshev_spline(n, k, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   x = x';
   % Caclulates data values at x 
   expx = exp(x);

   sp = spapi(optknt(x,k), x, expx);

   yplot = fnval(xplot,sp);

   % The function will calculate and display errors
   display_errors('Spline (Chebyshev)',n, yplot, expplot);
   
   %do_plotting_singular(yplot, x, expx);
end

% A function for interpolation using pchip 
% with n evenly spaced points
function [yplot] = evenly_pchip(n, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Caclulates data values at x 
   expx = exp(x);

   % Values using cubic spline
   yplot = interp1(x,expx, xplot,'pchip');

   % The function will calculate and display errors
   display_errors('Pchip (evenly)',n, yplot, expplot);
   
   %do_plotting_singular(yplot, x, expx);
end

% A function for interpolation using pchip 
% with n Chebyshev spaced points
function [yplot] = chebyshev_pchip(n, x)
   % Creates eveny spaced points
   xplot = linspace(0,2,1001);  
   % True values for the function
   expplot = exp(xplot);

   % Caclulates data values at x 
   expx = exp(x);

   % Values using cubic spline
   yplot = interp1(x,expx, xplot,'pchip');

   % The function will calculate and display errors
   display_errors('Pchip (Chebyshev)',n, yplot, expplot);
   
   %do_plotting_singular(yplot, x, expx);
end

function [] = doAllVandermonde()
   fprintf('=================================================\n')
   f = figure;
   f.Name = 'Interpolation, Vandermonde, Evenly Spaced Points';
   f.Position(1:4) = [200 200 1200 450];
   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Creates a mesh of x for the n points
      x = linspace(0,2,n)';

      yplot = evenly_spaced_vander(n, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end

   f = figure;
   f.Name = 'Interpolation, Vandermonde, Chebyshev Spaced Points';
   f.Position(1:4) = [200 200 1200 450];
   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Defines the Chebyshev data points on range [0-2]
      x = zeros(n,1);
      for k = 1:n
         x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
      end

      yplot = chebyshev_spaced_vander(n, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end

end

function [] = doAllBarycentric()
   fprintf('=================================================\n')   
   f = figure;
   f.Name = 'Interpolation, Barycentric, Evenly Spaced Points';
   f.Position(1:4) = [200 200 1200 450];
   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Creates a mesh of x for the n points
      x = linspace(0,2,n)';

      yplot = evenly_spaced_bary(n, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end

   f = figure;
   f.Name = 'Interpolation, Barycentric, Chebyshev Spaced Points';
   f.Position(1:4) = [200 200 1200 450];
   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Defines the Chebyshev data points on range [0-2]
      x = zeros(n,1);
      for k = 1:n
         x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
      end

      yplot = chebyshev_spaced_bary(n, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end
end

function [] = doAllSpline()
   fprintf('=================================================\n')   
   f = figure;
   f.Name = 'Interpolation, Spline, Evenly Spaced Points';
   f.Position(1:4) = [200 200 1200 450];
   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Creates a mesh of x for the n points
      x = linspace(0,2,n)';

      yplot = evenly_spline(n, 6, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end

   f = figure;
   f.Name = 'Interpolation, Spline, Chebyshev Spaced Points';
   f.Position(1:4) = [200 200 1200 450];
   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Defines the Chebyshev data points on range [0-2]
      x = zeros(n,1);
      for k = 1:n
         x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
      end
      yplot = chebyshev_spline(n, 6, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end
end

function [] = doAllPchip()
   fprintf('=================================================\n') 
   f = figure;
   f.Name = 'Interpolation, Pchip, Evenly Spaced Points';
   f.Position(1:4) = [200 200 1200 450];

   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Creates a mesh of x for the n points
      x = linspace(0,2,n)';
      yplot = evenly_pchip(n, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end

   f = figure;
   f.Name = 'Interpolation, Pchip, Chebyshev Spaced Points';
   f.Position(1:4) = [200 200 1200 450];

   i = 1;
   for n = [6,11,21,41,81,161,321,641]
      % Defines the Chebyshev data points on range [0-2]
      x = zeros(n,1);
      for k = 1:n
         x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
      end

      yplot = chebyshev_pchip(n, x);

      do_plotting_all(i , n , yplot, x);

      i = i + 1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doAllVandermonde()
doAllBarycentric()
doAllSpline()
doAllPchip()
