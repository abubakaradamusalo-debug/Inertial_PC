function [diffe,k]=Ex2TVC(x)
 
% x: initial vector on [0,1] grid
n = length(x); % number of discretization points

s = linspace(0, 1, n)';

xn = x;
xnm1 = xn;    

k=1;
tol=10^(-6);
N=3000;

diffe=5;
  tic;

% Open the file for writing
fileID = fopen('TVCEx2', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
    
    betan=1/(k+100);
    epsn=1/(k+1)^2;
    lam=0.1;
    
    alpn=0.5;
    if diffe==0
    alphan=alpn;
    else
    alphan=min(epsn/diffe,alpn);
    end
    
    wn=xn+alphan.*(xn-xnm1);
    
    % Define v(s) using interpolation
    
    v = @(s) interp1(linspace(0, 1, n), wn, s, 'linear', 'extrap');
    
    % Compute F(wn)
    
    F_wn = compute_Fv(v, 1); 

    pyn=wn-lam.*F_wn;

    yn=proj(pyn, 1);
     
    % Define new v(s) for yn
     
    v_new = @(s) interp1(linspace(0, 1, n), yn, s, 'linear', 'extrap');
    
    F_yn = compute_Fv(v_new, 1);
    
    
    dn=wn-yn-lam.*(F_wn-F_yn);
    thetan=0.8.*(norm(wn-yn).^2)/(norm(dn).^2);
     
    
    xnm1=xn;
    xn= betan*0.95*xnm1+(1-betan)*(wn-thetan.*dn);
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end
