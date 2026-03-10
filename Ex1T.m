function [err,k]=Ex1T(n)
   
xn= ones(n,1);    
    
j=1;
    
while(j <=n)
        
    xn(j)= (2*j+1)/(j^2+1);
        
    j=j+1;
    
end
     
ynm1=xn;
  

lamn=0.01;  

N=3000;
tol=10^(-6);

k=1;
err=1;
  tic;

% Open the file for writing
fileID = fopen('E1Tong', 'w');

while (k <= N && err > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, err);
     
    fprintf(' %d %d\n', k, err );  
       
     
    pyn=xn-lamn*F(ynm1);
    
    yn=proj_box(pyn,-2,2);   
     
    xnp1=yn+lamn*(F(ynm1)-F(yn)); 
     
    
    
    err= norm(xnp1-xn);
           
    xn=xnp1;
    ynm1=yn;
    k=k+1;
    
    
    
end
 
      toc; 
    



 end