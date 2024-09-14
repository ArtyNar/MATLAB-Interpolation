clear all
close all      
warning('off') 

% A helper function that estimates errors
 %   n:               number of points used
 %   yplot:           estimated values
 %   expplot:         true values
 %   errinf :          returns the infinity error found for future plotting
function [errinf] = display_errors(n,yplot, expplot)
    % Estimating errors
    errinf = norm((yplot-expplot),inf);     % estimate of infinity error 

    %fprintf('n = %3i |  inf error = %8.2e \n', n,errinf);
end

% A function that can be used to display a single graph
 %     inf_error:   infinity error found
 %     x:       set of x values at points used for estimation 
 %     y:       set of y values at points used for estimation
function [errinf] = interp_vander(n, x)
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
 
    % The function will calculate, display, and return errors
    errinf = display_errors(n, yplot, expplot);
 
    % Uses a function to simplify plotting
    % do_plotting_singular(yplot, x, expx); 
end

% A function for Polynomial Lagrange Interpolation using bary/baryweights routines 
 %   inf_error:   infinity error found
 %   n : number of points used to interpolate
 %   x : a set of points used to interpolate (evenly or chebyshev)
function [errinf] = interp_bary(n, x)
    % Creates eveny spaced points
    xplot = linspace(0,2,1001);  
    % True values for the function
    expplot = exp(xplot);
 
    % Caclulates data values at x 
    expx = exp(x);     
    
    vk = baryWeights(x); % Barycentric weights
    yplot = bary(xplot, expx, x, vk);
 
    % The function will calculate and display errors
    errinf = display_errors(n, yplot, expplot);
 
    %do_plotting_singular(yplot, x, expx);
end

% A function for spline interpolation 
 %   inf_error:   infinity error found
 %   n : number of points used to interpolate
 %   x : a set of points used to interpolate (evenly or chebyshev)
 %    k : order of spline
function [errinf] = interp_spline(n, k, x)
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
    errinf = display_errors(n, yplot, expplot);
    
    %do_plotting_singular(yplot, x, expx);
end

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


% A function for interpolation using pchip 
 %   inf_error:   infinity error found
 %   n : number of points used to interpolate
 %   x : a set of points used to interpolate (evenly or chebyshev)
function [errinf] = interp_pchip(n, x)
    % Creates eveny spaced points
    xplot = linspace(0,2,1001);  
    % True values for the function
    expplot = exp(xplot);
 
    % Caclulates data values at x 
    expx = exp(x);
 
    % Values using cubic spline
    yplot = interp1(x,expx, xplot,'pchip');
 
    % The function will calculate and display errors
    errinf = display_errors(n, yplot, expplot);
    
    %do_plotting_singular(yplot, x, expx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% A set of n points used in interolation
x_vals = [6,11,21,41,81,161,321,641];

% Pre-allocating for speed
errors_y_vander = zeros(1,8);
vand_time = zeros(1,8);
errors_y_bary = zeros(1,8);
bary_time = zeros(1,8);
errors_y_spline = zeros(1,8);
spline_time = zeros(1,8);
errors_y_pchip = zeros(1,8);
pchip_time = zeros(1,8);

% Used to index errors_y
i=1;

% Runs the algorithms and collects their errors
for n = [6,11,21,41,81,161,321,641]
    % Defines the Chebyshev data points on range [0-2]
    x = zeros(n,1);
    for k = 1:n
       x(k) = 1 - cos((k-1)*pi/(n-1)); % Chebyshev data points on range [0,2]
    end

    % Collects the errors into arrays and time it
    tic
        errors_y_vander(i) = interp_vander(n, x);
    vand_time(i) = toc;
    tic
        errors_y_bary(i) = interp_bary(n, x);
    bary_time(i) = toc;
    tic
        errors_y_spline(i) = interp_spline(n, 6, x);
    spline_time(i) = toc;
    tic
        errors_y_pchip(i) = interp_pchip(n, x);
    pchip_time(i) = toc;

    i = i + 1;
end




f = figure;
f.Name = 'Errors (Chebyshev)';
semilogy(x_vals, errors_y_vander,'-o',x_vals, errors_y_bary,'-o', x_vals, errors_y_spline,'-o', x_vals, errors_y_pchip,'-o', 'MarkerSize',3, 'LineWidth', 1)

% Adds legend to the plot
title('Errors (Chebyshev)')
legend('Vandermonde', 'Lagrange', 'Spline', 'Pchip');
xlabel(' n points');
ylabel('Error');

f = figure;
f.Name = 'Times (Chebyshev)';
plot(x_vals, vand_time,'-o', x_vals, bary_time,'-o', x_vals, spline_time,'-o', x_vals, pchip_time,'-o','MarkerSize',3, 'LineWidth', 1)

% Adds legend to the plot
title('Times (Chebyshev)')
legend('Vandermonde', 'Lagrange', 'Spline', 'Pchip');
xlabel('n points');
ylabel('Time (seconds)');

%%%%%
% Now, to repeat the same process fof the evenly-spaced points...
%%%%%

% Pre-allocating for speed
errors_y_vander = zeros(1,8);
vand_time = zeros(1,8);
errors_y_bary = zeros(1,3);
bary_time = zeros(1,3); % Note: ommiting when the values blow up
errors_y_spline = zeros(1,8);
spline_time = zeros(1,8);
errors_y_pchip = zeros(1,8);
pchip_time = zeros(1,8);

% Used to index errors_y
i=1;

% Runs the algorithms and collects their errors
for n = [6,11,21,41,81,161,321,641]
    % Defines the evenly data points on range [0-2]
    x = linspace(0,2,n)';  

    % Collects the errors into arrays and time it
    tic
        errors_y_vander(i) = interp_vander(n, x);
    vand_time(i) = toc;
    
    % Ommitting when the values blow up
    if ismember(n, [4 11 21 41])
        tic
            errors_y_bary(i) = interp_bary(n, x);
        bary_time(i) = toc;
    end
    tic
        errors_y_spline(i) = interp_spline(n, 6, x);
    spline_time(i) = toc;
    tic
        errors_y_pchip(i) = interp_pchip(n, x);
    pchip_time(i) = toc;

    i = i + 1;
end

f = figure;
f.Name = 'Errors (Evenly)';
semilogy(x_vals, errors_y_vander,'-o',[6 11 21 41], errors_y_bary,'-o', x_vals, errors_y_spline,'-o', x_vals, errors_y_pchip,'-o', 'MarkerSize',3, 'LineWidth', 1)

% Adds legend to the plot
title('Errors (Evenly)')
legend('Vandermonde', 'Lagrange', 'Spline', 'Pchip');
xlabel(' n points');
ylabel('Error');

f = figure;
f.Name = 'Times (Evenly)';
plot(x_vals, vand_time,'-o', [6 11 21 41], bary_time,'-o', x_vals, spline_time,'-o', x_vals, pchip_time,'-o','MarkerSize',3, 'LineWidth', 1)

% Adds legend to the plot
title('Times (Evenly)')
legend('Vandermonde', 'Lagrange', 'Spline', 'Pchip');
xlabel('n points');
ylabel('Time (seconds)');