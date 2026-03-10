function [diffe,k]=Ex2DCZR(x)
 
% x: initial vector on [0,1] grid
n = length(x); % number of discretization points

s = linspace(0, 1, n)';

xn = x;
xnm1 = xn;                   % Initialize xnm1



k=1;
tol=10^(-6);
N=3000;

diffe=5;
  tic;

% Open the file for writing
fileID = fopen('DCZREx2', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
    
    alpn=0.9*(k-1)/(k);  
    tau=0.9;
    gama=0.5;
     
    
    wn=xn+alpn.*(xn-xnm1);
    

    % Define v(s) using interpolation
    
    v = @(s) interp1(linspace(0, 1, n), wn, s, 'linear', 'extrap');
    
    % Compute F(wn)
    
    F_wn = compute_Fv(v, 1); 

    pyn=wn-tau.*F_wn;

    yn=proj(pyn, 1);
     
    % Define new v(s) for yn
     
    v_new = @(s) interp1(linspace(0, 1, n), yn, s, 'linear', 'extrap');
    
    F_yn = compute_Fv(v_new, 1);
    
    dn=wn-yn-tau.*(F_wn-F_yn);
    
    con1=norm(dn);
    if con1==0    
        betan=0;
    else
        betan=inner_L2(wn-yn,dn,s)/(con1.^2);
    end
     
    
    xnm1=xn;
    xn= wn-gama*betan.*dn;
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end
