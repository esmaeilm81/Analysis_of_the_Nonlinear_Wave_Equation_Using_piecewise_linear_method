clc;
clear;
close all;

N=1001;
h=10/(N-1);
u_max=1;
xi=-4.995:0.01:4.995;
x=linspace(-5,5,1000);
a=length(xi);
u0=zeros(1,a);
u0(x>-1 & x<0)=-1;
u0(x>0 & x<1)=1;
u0(2,:)=u0(1,:);
u0(3,:)=u0(1,:);
t=0.5*h/u_max;
nt=floor(2/t);

u_p=u0;
f=@(u) 0.5*u^2
Mind=@(r) max(0,min(r,1));
sbee=@(r) max(max(0,min(2*r,1)),min(r,2));
non=@(r) 1;
Lim={non Mind sbee};
for j=1:3
    for n=1:nt
        u3(j,:)=u_p(j,:);
        for i=3:a-2
            uL_1=U(u_p(j,i-1),x(i-1)+0.5*h,x(i-1),Lim{j},r(u_p(j,i-2),u_p(j,i-1),u_p(j,i)),s(u_p(j,i),u_p(j,i-1),h));
            uR_1=U(u_p(j,i),x(i-1)+0.5*h,x(i),Lim{j},r(u_p(j,i-1),u_p(j,i),u_p(j,i+1)),s(u_p(j,i+1),u_p(j,i),h));
            uL_2=U(u_p(j,i),x(i)+0.5*h,x(i),Lim{j},r(u_p(j,i-1),u_p(j,i),u_p(j,i+1)),s(u_p(j,i+1),u_p(j,i),h));;
            uR_2=U(u_p(j,i+1),x(i)+0.5*h,x(i+1),Lim{j},r(u_p(j,i),u_p(j,i+1),u_p(j,i+2)),s(u_p(j,i+2),u_p(j,i+1),h));
            flux_left=R_flux(uL_1, uR_1, f);
            flux_right=R_flux(uL_2, uR_2, f);
            u3(j,i)=u_p(j,i)-(t/h)*(flux_right - flux_left);
        end
        u_p(j,:)=u3(j,:);
    end
end
% Plot the result
plot(x, u3, 'LineWidth', 2);
title('Solution at t = 2 s using Roe''s Method');
xlabel('x');
ylabel('u(x, t)');


function F=R_flux(uL,uR,f)
    F=(uL>uR)*(max(f(uL),f(uR)))+...
    (uL<=uR)*(min(f(uL),f(uR)));
end
function a=r(uL,uM,uR)
    a=(uM-uL)/(uR-uM);
end
function a=U(ui,x ,xi,phi,r,s)
    a=ui+phi(r)*s*(x-xi)/2;
end
function a=s(uR,uL,h)
    a=(uR-uL)/h;
end