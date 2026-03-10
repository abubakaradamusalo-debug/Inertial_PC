function [diffe,k]=Ex2CTC(x)
 

% x: initial vector on [0,1] grid
n = length(x); % number of discretization points

s = linspace(0, 1, n)';

xn = x;
xnm1 = xn;                   % Initialize xnm1
     
alpha=0.9;  
lamn=0.17;
gama=1.9;

k=1;
tol=10^(-6);
N=3000;

diffe=5;
  tic;

% Open the file for writing
fileID = fopen('CTCEx2', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
     
    taun=1/(k+1)^3; 
    thetan=1/(k+1000)^0.5;
    betan=0.9*(1-thetan);



    if norm(xn-xnm1)==0
        alphan=alpha;
    else
        alphan= min(taun/norm(xn-xnm1), alpha);
    end
     
    
    wn=xn+alphan*(xn-xnm1);
    

    % Define v(s) using interpolation
    
    v = @(s) interp1(linspace(0, 1, n), wn, s, 'linear', 'extrap');
    
    % Compute F(wn)
    
    F_wn = compute_Fv(v, 1); 

    pyn=wn-lamn.*F_wn;
    yn=proj(pyn, 1);
    
    % Define new v(s) for yn
     
    v_new = @(s) interp1(linspace(0, 1, n), yn, s, 'linear', 'extrap');
    
    F_yn = compute_Fv(v_new, 1);
    
    dn=wn-yn-lamn*(F_wn-F_yn);
    
    con1=norm(dn);
    if con1==0    
        etan=0;
    else
        etan=inner_L2(wn-yn,dn,s)/(con1.^2);
    end
     

    zn=wn-gama*lamn*etan*dn;  

    xnp1=(1-thetan-betan)*wn+betan*zn;


    xnm1=xn;
    xn=xnp1; 
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end
