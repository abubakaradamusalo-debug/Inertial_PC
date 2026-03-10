function [Tcpu, Itr, NormF, xn] = IDFPM(fun, method, Switch, model, x0, para)
    format long;
    tic;
    epsilon = 1e-6;
    k_max = para.Itr_max;
    gamma = para.gamma;
    sigma = para.sigma;
    tau = para.tau;
    rho = para.rho;
    err = 1;
    xnm2 = x0;
    xnm1=xnm2;
    xn = xnm1; 
    ynm1=x0;

    lamn = 0.25; 
                
    thetanm1=1;
    k = 1;

    while (k <= k_max && err > epsilon)
        switch method 
            case 'PASYAlg3p3'
 

                phi0=0.01;
                phi1=0.35;
 
                thetan=-0.25*(k+1)/(k);
                taun=1/(k+1)^3; 
                gama=0.01; 

                wn=(xn)+thetan*((xn)-(xnm1)); 

                Fynm1 = feval(fun, ynm1);
                pyn = wn-lamn*Fynm1;
                yn = soft_threshold(pyn, 0.25*lamn);
                Fyn = feval(fun, yn);                 
            
                un=yn-lamn*(Fynm1-Fyn);
            
                dn=wn-yn-lamn*(Fynm1-Fyn);
                             
                if dn==0
                    betan=0;
                else
                    betan=dot(un-yn,dn)/norm(dn)^2;
                end
                 
                
                xnp1=un-gama*betan*dn; 
                
                con= norm(Fynm1-Fyn);
                
                if con>(phi0/lamn)*norm(ynm1-yn)
                    lamnp1=phi1*norm(ynm1-yn)/con;
                else
                    lamnp1=(1+taun)*lamn;
                end
                 

                ynm1=yn;
                 

            case 'CTCAlg3p1'

                
                alpha=0.9;  
                gama=1.9;

                taun=1/(k+1)^3; 
                thetan=1/(k+1000);
                betan=0.9*(1-thetan);
            
                if norm(xn-xnm1)==0
                    alphan=alpha;
                else
                    alphan= min(taun/norm(xn-xnm1), alpha);
                end
                 

                wn=(xn)+alphan*((xn)-(xnm1)); 
                
                
                Fwn = feval(fun, wn);
                pyn = wn-lamn*Fwn;
                yn = soft_threshold(pyn, 0.25*lamn);
                Fyn = feval(fun, yn);
                 
                            
                dn=wn-yn-lamn*(Fwn-Fyn);
                
                con1=norm(dn);
                if con1==0    
                    etan=0;
                else
                    etan=dot(wn-yn,dn)/(con1.^2);
                end
                 
            
                zn=wn-gama*lamn*etan*dn;  
            
                xnp1=(1-thetan-betan)*wn+betan*zn;

 
                lamnp1=lamn;



            case 'DCZRAlg3p1'
                
            
                alpn=0.9*(k-1)/(k);   
                gama=0.1;
                 
                
                wn=xn+alpn.*(xn-xnm1); 


                Fwn = feval(fun, wn);
                pyn = wn-lamn*Fwn;
                yn = soft_threshold(pyn, 0.25*lamn);
                Fyn = feval(fun, yn);
                             
                dn=wn-yn-tau.*(Fwn-Fyn);
                
                con1=norm(dn);
                if con1==0    
                    betan=0;
                else
                    betan=dot(wn-yn,dn)/(con1.^2);
                end
                  
                
                xnp1=wn-gama*betan.*dn;                 
                  
                lamnp1=lamn;


                 
            case 'MAlg1'
                             
                blam=0.001; 
                phi=1.7;  
                 
                rho=1/phi+1/phi^2;
                
                e1=rho*lamn;

                Fxnm1 = feval(fun, xnm1);
                Fxn = feval(fun, xn);

                e2=(phi*thetanm1/(4*lamn))*norm(xnm1-xn)/norm(Fxnm1-Fxn);
                
                E=[e1;e2;blam];
                lamnp1= min(E);
                
                bzk=(1/phi)*((phi-1)*xn+ynm1); 
                
                pxnp1=bzk-lamnp1*Fxn;
                
                xnp1=soft_threshold(pxnp1,0.25*lamnp1);  
                
                thetak=(lamnp1/lamn)*phi;
                
                thetanm1=thetak;
                   
                   
            case 'TQYAlg3p1'
          
                mu=0.9;
                tau=0.1;
                theta=1.5;
                epsn=1/(k+1)^2;
                zetan=1/(k+1)^3;
                sigman=1/(k+1000);
                psin=0.35*(1-sigman);
            
                if norm(xn-xnm1)==0
                    taun=tau;
                else
                    taun= min(epsn/norm(xn-xnm1), tau);
                end
                  
                wn=(xn)+taun*((xn)-(xnm1)); 

                Fwn = feval(fun, wn);
                pyn = wn-lamn*Fwn;
                yn = soft_threshold(pyn, 0.25*lamn);
                Fyn = feval(fun, yn);
                            
             
                cn=wn-yn-lamn*(Fwn-Fyn);
                
                con1=norm(cn);
                if con1==0    
                    chin=0;
                else
                    chin=dot(wn-yn,cn)/(con1.^2);
                end
             
                pzn=wn-theta*lamn*chin*Fyn;
                nn=wn-lamn*Fwn-yn;
                zn=proj_halfspace(pzn, nn, yn);
            
                con=norm(Fwn-Fyn);
                
                if con==0
                    lamnp1=lamn+zetan;
                else
                    lamnp1=min(mu*norm(wn-yn)/con, lamn+zetan);
                end
            
                xnp1=(1-sigman-psin)*wn+psin*zn;

                
                

            case 'TVCAlg1'
                
            
                betan=1/(k+100);
                epsn=1/(k+1)^2; 
                
                alpn=0.5;
                if norm(xn-xnm1)==0
                alphan=alpn;
                else
                alphan=min(epsn/norm(xn-xnm1),alpn);
                end
                 
                
                wn=(xn)+alphan*((xn)-(xnm1)); 
                
                Fwn = feval(fun, wn);
                pyn = wn-lamn*Fwn;
                yn = soft_threshold(pyn, 0.25*lamn);
                Fyn = feval(fun, yn); 

                dn=wn-yn-lamn*(Fwn-Fyn);
                thetan=0.8.*(norm(wn-yn).^2)/(norm(dn).^2);
                
                xnp1= betan*0.95*xnm1+(1-betan)*(wn-thetan.*dn); 

                lamnp1 = lamn; 

            case 'TAlg3p1'
 

                Fynm1 = feval(fun, ynm1);
                pyn = xn-lamn*Fynm1;
                yn = soft_threshold(pyn, 0.25*lamn);
                Fyn = feval(fun, yn);                 
              

                xnp1=yn+lamn*(Fynm1-Fyn); 
    
                lamnp1 =  lamn; 


            otherwise
                disp('Input error! Please check the input method');
        end
        xnm1 = xn;
        xn = xnp1;
        lamn = lamnp1;
        err = norm(xn - xnm1);
        k = k + 1;
    end
    Itr = k;
    Tcpu = toc;
    NormF = err;
end
