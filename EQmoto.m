function [theta,a,b,g,d,tempo,thp]=EQmoto(t,a,b,g,d,a0,b0,c0,d0,e0,l0,DL0,L0,m1,m2,m3,J0,J3,k,r,C,tp0,h,Tfinal)

eps=0.003;

%calculating some parameters needed in the motion equations
JA = @(t,a,b,g,d) J0 + (m1*a0^2)+(m2*a0^2)*((sin(a-t))/(sin(a-b)))^2 +
(m3*(a0*l0/e0)^2+J3*(a0/e0)^2)*((sin(a-t))/(sin(a-b)))^2*((sin(g-b))/(sin(g-d)))^2;
dJdT= @(t,a,b,g,d) m2*a0*((sin(2*a-2*t)-sin(2*a-2*b))*a0*sin(b-t)/(((sin(ab))^4)*b0*sin(a-b))
- sin(2*a-2*t)/((sin(a-b))^2) + ((sin(a-t))^2)*cos(a-b)*a0*(sin(at))/((sin(a-b))*c0*sin(a-b)))
+ (m3*(a0*l0/e0)^2 + J3*(a0/e0)^2)*((sin(2*a-2*t)-sin(2*a-
2*b))*a0*sin(b-t)/(((sin(a-b))^4)*b0*sin(a-b)) - sin(2*a-2*t)/((sin(a-b))^2) + ((sin(at))^2)*cos(a-b)*a0*(sin(a-t))/((sin(a-b))*c0*sin(a-b)))*((sin(g-b))^2)/((sin(g-d))^2)
+
((sin(a-t))^2)/((sin(a-b))^2)*((sin(2*g-2*b)-sin(2*g-2*d))*a0*sin(a-t)*sin(d-b)/(((sin(gd))^4)*d0*sin(a-b)*sin(g-d))
- sin(2*g-2*b)/((sin(g-d))^2)*(a0/c0*sin(a-t)/sin(a-b)) +
((sin(g-b))^2)*cos(g-d)*a0*(sin(a-t))*sin(g-b)/((sin(g-d))*e0*sin(a-b)*sin(g-d)));
dvdthetaELA= @(t,a,b,g,d,tempo) k(tempo)*(-L0+DL0 + (a0*cos(t)+ b0*cos(a)+ d0*cos(g)))*(-
a0*sin(t)-b0*sin(a)*a0*sin(b-t)/(b0*sin(a-b)) - sin(g)*a0*sin(a-t)*sin(d-b)/(sin(ab)*sin(g-d)));
dvdthetaGRA= @(t,a,b,g,d) 9.8*( (m1+m2+m3)*a0*cos(t) + (m2+m3)*a0*(sin(at)*cos(b))/(sin(a-b))
+ m3*l0*a0*(sin(a-t))/(sin(a-b))*(sin(g-b))/(sin(g-d))/e0*cos(d));
dRdthetap= @(t,a,b,g,d,tempo) r(tempo)*(((l0*a0/e0)^2)*(((sin(a-t))/(sin(ab)))^2)*(((sin(g-b))/(sin(g-d)))^2))*sin(d)^2;

%calculating Jacobian elements
Jalpha= @(t,a,b,g,d) a0/b0*sin(b-t)/sin(a-b);
Jbeta= @(t,a,b,g,d) a0/c0*sin(a-t)/sin(a-b);
Jgamma= @(t,a,b,g,d) a0/d0*sin(a-t)/sin(a-b)*(sin(d-b))/(sin(g-d));
Jdelta= @(t,a,b,g,d) a0/e0*sin(a-t)/sin(a-b)*(sin(g-b))/(sin(g-d));

%Euler's method
tempo=0:h:Tfinal;
theta(1)=t;
tp=tp0;
thp(1)=tp0;

for i=1:length(tempo)-1
 acc=-(0.5*dJdT(theta(i),a(i),b(i),g(i),d(i))*tp^2+(dvdthetaELA(theta(i),a(i),b(i),g(i),d(i),tempo(i)) + dvdthetaGRA(theta(i),a(i),b(i),g(i),d(i))+
      dRdthetap(theta(i),a(i),b(i),g(i),d(i),tempo(i))*tp -C))/JA(theta(i),a(i),b(i),g(i),d(i));
 tp1=tp+acc*h;
 theta(i+1)=theta(i)+tp*h+0.5*acc*h^2;
 a(i+1)=Jalpha(theta(i),a(i),b(i),g(i),d(i))*(theta(i+1)-theta(i))+ a(i);
 b(i+1)=Jbeta(theta(i),a(i),b(i),g(i),d(i))*(theta(i+1)-theta(i))+ b(i);
 g(i+1)=Jgamma(theta(i),a(i),b(i),g(i),d(i))*(theta(i+1)-theta(i))+g(i);
 d(i+1)=Jdelta(theta(i),a(i),b(i),g(i),d(i))*(theta(i+1)-theta(i)) + d(i);
 tp=tp1;
 thp(i+1)=tp;
end
