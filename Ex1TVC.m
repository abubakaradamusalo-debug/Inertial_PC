function [diffe,k]=Ex1TVC(n)


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
fileID = fopen('TVCEx1', 'w');

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
    
    pyn=wn-lam.*F(wn);
    yn=proj_box(pyn, -2,2);
     
    
    dn=wn-yn-lam.*(F(wn)-F(yn));
    thetan=0.8.*(norm(wn-yn).^2)/(norm(dn).^2);
     
    
    xnm1=xn;
    xn= betan*0.95*xnm1+(1-betan)*(wn-thetan.*dn);
    
    diffe= norm(xn-xnm1);
    
    k=k+1;
    
    
    
end
 
      toc; 
    



 end
