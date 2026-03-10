function [diffe,k]=Ex1Alg3p3(n)
 
xnm1= ones(n,1);    
    
j=1;
    
while(j <=n)
        
    xnm1(j)= (2*j+1)/(j^2+1);
        
    j=j+1;
    
end
     
    
xn=xnm1; 
ynm1=xnm1;


lamn=0.01;

phi0=0.01;
phi1=0.35;
 

 
k=1;
tol=10^(-6);
N=3000;
 

diffe=5;
  tic;


% Open the file for writing
fileID = fopen('Ex1Alg3p3', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    fprintf(' %d %d\n', k, diffe ); 
    
    thetan=-0.25*(k+1)/(k);
    taun=1/(k+1)^3; 
    gama=0.01; 

    wn=(xn)+thetan*((xn)-(xnm1)); 
    pyn=(wn)-lamn*F(ynm1);
    yn=proj_box(pyn, -2,2);
    
    un=yn-lamn*(F(ynm1)-F(yn));

    dn=wn-yn-lamn*(F(ynm1)-F(yn));

    if dn==0
        betan=0;
    else
        betan=dot(un-yn,dn)/norm(dn)^2;
    end
     
    
    xnp1=un-gama*betan*dn; 
    
    con= norm(F(ynm1)-F(yn));
    
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