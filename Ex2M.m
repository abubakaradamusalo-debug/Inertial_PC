function [diffe,k]=Ex2M(x)

 
% x: initial vector on [0,1] grid
n = length(x); % number of discretization points

s = linspace(0, 1, n)';

zkm1 =  x;
        
zk=zkm1;    
bzkm1=zk;
 
 

thetakm1=1;
  

 
lamkm1=9e-1;
blam=0.9; 
k=1;
tol=10^(-6);
N=3000;


diffe=5;
  tic;

% Open the file for writing

fileID = fopen('MEx4Alg1', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    fprintf(' %d %d\n', k, diffe );  
    phi=1.7; 
     
    rho=1/phi+1/phi^2;
    
    e1=rho*lamkm1;


    % Define v(s) using interpolation
    v = @(s) interp1(linspace(0, 1, n), zkm1, s, 'linear', 'extrap');
    
        % Compute F(zkm1)
    
    F_zkm1 = compute_Fv(v, 1); 


    % Define v(s) using interpolation
    v_zk = @(s) interp1(linspace(0, 1, n), zk, s, 'linear', 'extrap');
    
        % Compute F(zkm1)
    
    F_zk = compute_Fv(v_zk, 1); 


    e2=(phi*thetakm1/(4*lamkm1))*norm(zkm1-zk)/norm(F_zkm1-F_zk);
    
    E=[e1;e2;blam];
    lamk= min(E);
    
    bzk=(1/phi)*((phi-1)*zk+bzkm1); 
    
    pznp1=bzk-lamk*F_zk;
    
    znp1= proj(pznp1,1);  
    
    thetak=(lamk/lamkm1)*phi;
    
    
    diffe= norm(znp1-zk);
    
    zkm1=zk;  
    zk=znp1;    
    lamkm1=lamk;
    thetakm1=thetak;
    
    k=k+1;
    
    
end
      toc; 
    
 
 end
