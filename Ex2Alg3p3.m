function [diffe,k]=Ex2Alg3p3(x)

% x: initial vector on [0,1] grid
n = length(x); % number of discretization points


s = linspace(0, 1, n)';

xn = x;
xnm1 = xn;                   % Initialize xnm1
     
ynm1=xnm1;


lamn=0.1;

phi0=0.01;
phi1=0.35;
 

 
k=1;
tol=10^(-6);
N=3000;
 

diffe=5;
  tic;


% Open the file for writing
fileID = fopen('Ex2Alg3p3', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    fprintf(' %d %d\n', k, diffe ); 
    
    thetan=-1/(k+1000);
    taun=1/(k+1)^3; 
    gama=0.01; 

    wn=(xn)+thetan*((xn)-(xnm1)); 

    % Define v(s) using interpolation
    
    v = @(s) interp1(linspace(0, 1, n), ynm1, s, 'linear', 'extrap');
    
    % Compute F(wn)
    
    F_ynm1 = compute_Fv(v, 1); 

    pyn=(wn)-lamn*F_ynm1;

    yn=proj(pyn, 1);
    
    % Define new v(s) for yn
     
    v_new = @(s) interp1(linspace(0, 1, n), yn, s, 'linear', 'extrap');
    
    F_yn = compute_Fv(v_new, 1);

    un=yn-lamn*(F_ynm1-F_yn);

    dn=wn-yn-lamn*(F_ynm1-F_yn);

    if dn==0
        betan=0;
    else
        betan=inner_L2(un-yn,dn,s)/norm(dn)^2;
    end
     
    
    xnp1=un-gama*betan*dn; 
    
    con= norm(F_ynm1-F_yn);
    
    if con>(phi0/lamn)*norm(ynm1-yn)
        lamnp1=phi1*norm(ynm1-yn)/con;
    else
        lamnp1=(1+taun)*lamn;
    end
    
    lamn=lamnp1;  
    ynm1=yn;
    xnm1=xn;
    xn=xnp1;
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end