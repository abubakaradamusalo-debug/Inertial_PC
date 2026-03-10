function pj=proj(x,rad) 

xx0=  x;
cc=norm(xx0);
if cc<=rad
    pj=  xx0;
else
    pj=  (rad/cc)*xx0;
end
end
