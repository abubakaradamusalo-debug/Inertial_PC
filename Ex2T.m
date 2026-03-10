function [err,k]=Ex2T(x)
 
% x: initial vector on [0,1] grid
n = length(x); % number of discretization points

s = linspace(0, 1, n)';

xn =   x;

ynm1=xn;
     

lamn=0.1;  

N=5000;
tol=10^(-6);

k=1;
err=1;
  tic;

% Open the file for writing
fileID = fopen('E2T', 'w');

while (k <= N && err > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, err);
     
    fprintf(' %d %d\n', k, err );  
       

    % Define v(s) using interpolation
    v = @(s) interp1(linspace(0, 1, n), ynm1, s, 'linear', 'extrap');
    
        % Compute F(xn)
    
    F_ynm1 = compute_Fv(v, 1); 

     
    pyn=xn-lamn*F_ynm1;
    
    yn=proj(pyn, 1);

    % Define new v(s) for yn
     
    v_new = @(s) interp1(linspace(0, 1, n), yn, s, 'linear', 'extrap');
    
    F_yn = compute_Fv(v_new, 1);
    
     
    xnp1=yn+lamn*(F_ynm1-F_yn); 
     
    
    
    err= norm(xnp1-xn);
           
    xn=xnp1;
    k=k+1;
    
    
    
end
 
      toc; 
    




 end