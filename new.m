%clear all
%close all
 
warning('off')

% Creates eveny spaced points for error sampling
xplot = linspace(0,2,1001);  
% Calculates the true values for the function for error sampling
expplot = exp(xplot);

% This is a helper function that estimates errors
function [] = display_errors(n,yplot, expplot)
    % Estimating errors
    sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
    errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
    err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm

    fprintf('Vander (evenly-spaced) | n = %3i |  inf error = %8.2e | 2 norm error = %8.2e \n', n,errinf,err2);
end

% A function for Polynomial Interpolation using Vandermonde matrix 
% with n evenly spacedpoints
function [] = vander_interpolation(xplot, expplot)
    % Setting up a figure for future graphing
    f = figure;
    f.Name = 'Interpolation';
    f.Position(1:4) = [200 200 900 600];
    
    % Plots x mesh, true values of y
    subplot(3,3,1);
    plot(xplot, expplot, 'r');
 
    % Adds legend to the plot
    title('True graph')
    xlabel(' x ');
    ylabel('exp(x)');
 
    i = 2;
 
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

       % The function will calculate and display errors
       display_errors(n, yplot, expplot)

       % Plots x mesh, approximated values of y, and approximation points
       subplot(3,3,i);
       plot(xplot, yplot, x, expx, '*');
       
       % Adds legend to the plot
       title(strcat('Interpolated with n = ', num2str(n)))
       legend('Interpolated graph', 'Interpolation points');
       xlabel(' x ');
       ylabel('exp(x) - Approximated');
 
       i = i + 1;
    end
end

vander_interpolation(xplot, expplot)