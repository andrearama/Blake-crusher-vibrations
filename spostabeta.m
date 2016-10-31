function [t,a,b,g,d]= spostabeta(t,a,b,g,d,a0,b0,c0,d0,e0,l0,thetaspostato,h)

if thetaspostato>=t
x=t:h:thetaspostato;
else
x=t:-h:thetaspostato;
end

Jalpha= @(t,a,b,g,d) a0/b0*sin(b-t)/sin(a-b);
Jbeta= @(t,a,b,g,d) a0/c0*sin(a-t)/sin(a-b);
Jgamma= @(t,a,b,g,d) a0/d0*sin(a-t)/sin(a-b)*(sin(d-b))/(sin(g-d));
Jdelta= @(t,a,b,g,d) a0/e0*sin(a-t)/sin(a-b)*(sin(g-b))/(sin(g-d));

for i=2:length(x)
 dtheta=x(i)-t;
 a=Jalpha(t,a,b,g,d)*dtheta+ a;
 b=Jbeta(t,a,b,g,d)*dtheta+ b;
 g=Jgamma(t,a,b,g,d)*dtheta+ g;
 d=Jdelta(t,a,b,g,d)*dtheta+ d;
 t=x(i);
end
