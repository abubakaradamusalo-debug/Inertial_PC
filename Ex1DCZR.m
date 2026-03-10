function [diffe,k]=Ex1DCZR(n)
 
xnm1= ones(n,1);    
    
j=1;
    
while(j <=n)
        
    xnm1(j)= (2*j+1)/(j^2+1);
        
    j=j+1;
    
end
     
    
xn=xnm1; 
 



k=1;
tol=10^(-6);
N=3000;

diffe=5;
  tic;

% Open the file for writing
fileID = fopen('DCZREx1', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
    
    alpn=0.9*(k-1)/(k);  
    tau=0.1;
    gama=0.1;
     
    
    wn=xn+alpn.*(xn-xnm1);
    
    pyn=wn-tau.*F(wn);
    yn=proj_box(pyn, -2,2);
     
    
    dn=wn-yn-tau.*(F(wn)-F(yn));
    
    con1=norm(dn);
    if con1==0    
        betan=0;
    else
        betan=dot(wn-yn,dn)/(con1.^2);
    end
     
    
    xnm1=xn;
    xn= wn-gama*betan.*dn;
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end
