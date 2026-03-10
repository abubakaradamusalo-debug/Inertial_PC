function [diffe,k]=Ex1M(n)

 

zkm1= ones(n,1);    
    
j=1;
    
while(j <=n)
        
    zkm1(j)= (2*j+1)/(j^2+1);
        
    j=j+1;
    
end  
zk=zkm1;
bzkm1=zk;
 
thetakm1=1;
 
 
 
lamkm1=0.1;
blam=0.001; 
k=1;
tol=10^(-6);
N=3000;


diffe=5;
  tic;

% Open the file for writing

fileID = fopen('MEx1Alg1', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    fprintf(' %d %d\n', k, diffe );  
    phi=1.7; 
     
    rho=1/phi+1/phi^2;
    
    e1=rho*lamkm1;
    e2=(phi*thetakm1/(4*lamkm1))*norm(zkm1-zk)/norm(F(zkm1)-F(zk));
    
    E=[e1;e2;blam];
    lamk= min(E);
    
    bzk=(1/phi)*((phi-1)*zk+bzkm1); 
    
    pznp1=bzk-lamk*F(zk);
    
    znp1= proj_box(pznp1,-2,2);  
    
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
