run Dati
k=@(t) k0;
r=@(t) r0;

%Stable equilibrium positions calculus
x=linspace(0,2*pi,720);
theta01=theta0;
alpha01=alpha0;
beta01=beta0;
gamma01=gamma0;
delta01=delta0;

for i=1:length(x)
thetaspostato=x(i);
[theta0,alpha0,beta0,gamma0,delta0] =spostabeta(theta01,alpha01,beta01,gamma01,delta01,a0,b0,c0,d0,e0,l0,thetaspostato,0.001);
DL= @(t,a,b,g,d) (C - 9.8*( m1*a0*cos(t) + m2*a0*(sin(a-t)*cos(b))/(sin(a-b)) +m3*l0*a0*(sin(a-t))/(sin(a-b))*(sin(g-b))/(sin(g-d))/e0*cos(d)))/(-a0*sin(t)-
    b0*sin(a)*a0*sin(b-t)/(b0*sin(a-b)) - sin(g)*a0*sin(a-t)*sin(d-b)/(sin(a-b)*sin(g-d)))/k0- (a0*cos(t)+ b0*cos(a)+ d0*cos(g)) + L0;
DL0(i)=DL(theta0,alpha0,beta0,gamma0,delta0);

end

plot(x,DL0)
