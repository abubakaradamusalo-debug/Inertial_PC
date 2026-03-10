function [diffe,k]=Ex1CTC(n)
 
xnm1= ones(n,1);    
    
j=1;
    
while(j <=n)
        
    xnm1(j)= (2*j+1)/(j^2+1);
        
    j=j+1;
    
end
     
    
xn=xnm1; 

alpha=0.9;  
lamn=0.17;
gama=1.9;

k=1;
tol=10^(-6);
N=3000;

diffe=5;
  tic;

% Open the file for writing
fileID = fopen('CTCEx1', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
     
    taun=1/(k+1)^3; 
    thetan=1/(k+1000);
    betan=0.9*(1-thetan);



    if norm(xn-xnm1)==0
        alphan=alpha;
    else
        alphan= min(taun/norm(xn-xnm1), alpha);
    end
     
    
    wn=xn+alphan*(xn-xnm1);
    
    pyn=wn-lamn.*F(wn);
    yn=proj_box(pyn, -2,2);
    
    
    dn=wn-yn-lamn*(F(wn)-F(yn));
    
    con1=norm(dn);
    if con1==0    
        etan=0;
    else
        etan=dot(wn-yn,dn)/(con1.^2);
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
