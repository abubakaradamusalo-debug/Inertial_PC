function [diffe,k]=Ex2TQY(x)
 
% x: initial vector on [0,1] grid
n = length(x); % number of discretization points

s = linspace(0, 1, n)';

xn = x;
xnm1 = xn;     

tau=0.1; 
mu=0.9;
nun=0.1;
theta=1.5;

k=1;
tol=10^(-6);
N=3000;

diffe=5;
  tic;

% Open the file for writing
fileID = fopen('TQYEx2', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
     
    epsn=1/(k+1)^2;
    zetan=1/(k+1)^3;
    sigman=1/(k+1);
    psin=0.9*(1-sigman);



    if norm(xn-xnm1)==0
        taun=tau;
    else
        taun= min(epsn/norm(xn-xnm1), tau);
    end
     
    
    un=xn+taun*(xn-xnm1);
    

    % Define v(s) using interpolation
    
    v = @(s) interp1(linspace(0, 1, n), un, s, 'linear', 'extrap');
    
    % Compute F(un)
    
    F_un = compute_Fv(v, 1); 

    pyn=un-nun.*F_un;
    yn=proj(pyn, 1);

    % Define new v(s) for yn
     
    v_new = @(s) interp1(linspace(0, 1, n), yn, s, 'linear', 'extrap');
    
    F_yn = compute_Fv(v_new, 1);
        
    cn=un-yn-nun*(F_un-F_yn);
    
    con1=norm(cn);
    if con1==0    
        chin=0;
    else
        chin=dot(un-yn,cn)/(con1.^2);
    end
     

    pzn=un-theta*nun*chin*F_yn;
    nn=un-nun*F_un-yn;
    zn=proj_halfspace(pzn, nn, yn);

    con=norm(F_un-F_yn);
    
    if con==0
        nunp1=nun+zetan;
    else
        nunp1=min(mu*norm(un-yn)/con, nun+zetan);
    end

    xnp1=(1-sigman-psin)*un+psin*zn;


    xnm1=xn;
    xn=xnp1;
    nun=nunp1;
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end
