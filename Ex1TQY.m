function [diffe,k]=Ex1TQY(n)
 
xnm1= ones(n,1);    
    
j=1;
    
while(j <=n)
        
    xnm1(j)= (2*j+1)/(j^2+1);
        
    j=j+1;
    
end
     
    
xn=xnm1; 

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
fileID = fopen('TQYEx1', 'w');

while (k <= N && diffe > tol)
    % Write output to file
    fprintf(fileID, '%d %d\n', k, diffe);
    
    
    fprintf(' %d %d \n', k, diffe);  
     
    epsn=1/(k+1)^2;
    zetan=1/(k+1)^3;
    sigman=1/(k+1000);
    psin=0.35*(1-sigman);



    if norm(xn-xnm1)==0
        taun=tau;
    else
        taun= min(epsn/norm(xn-xnm1), tau);
    end
     
    
    un=xn+taun*(xn-xnm1);
    
    pyn=un-nun.*F(un);
    yn=proj_box(pyn, -2,2);
    
    
    cn=un-yn-nun*(F(un)-F(yn));
    
    con1=norm(cn);
    if con1==0    
        chin=0;
    else
        chin=dot(un-yn,cn)/(con1.^2);
    end
     

    pzn=un-theta*nun*chin*F(yn);
    nn=un-nun*F(un)-yn;
    zn=proj_halfspace(pzn, nn, yn);

    con=norm(F(un)-F(yn));
    
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
