% This is a helper function that estimates errors
function [] = display_errors(n,yplot, expplot)
    % Estimating errors
    sqrth = 1.0/sqrt(n); % error factor for spatial l2 norm
    errinf = norm((yplot-expplot),inf);     % estimate of infinity error 
    err2   = sqrth*norm((yplot-expplot),2); % estimate of 2 norm

    fprintf('Vander (evenly-spaced) | n = %3i |  inf error = %8.2e | 2 norm error = %8.2e \n', n,errinf,err2);
end