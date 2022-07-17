function graphene_plot_band_structure_qc_debug(V,Vg,Vi,Q_top,Q_bottom,flag)

if nargin < 6
    flag = 0;
end

er = 3;
vf = 1e6;
e_el = 1.602e-19;
Vf = 1e6;
hbar = 1.055e-34;
Ef_top =- sqrt(4*pi*abs(Q_top)/e_el).*sign(Q_top)*hbar*Vf/(2*e_el);
Ef_bott = -sqrt(4*pi*abs(Q_bottom)/e_el).*sign(Q_bottom)*hbar*Vf/(2*e_el);

ktop = e_el*Ef_top/(hbar*Vf);
kbott= e_el*Ef_bott/(hbar*Vf);
k_i = abs(e_el*Vi/(hbar*Vf))
%x= linspace(-1.5*max(abs(ktop+k_i),abs(kbott+k_i)),1.5*max(abs(ktop+k_i),abs(kbott+k_i)),100);
x= linspace(-1.5*e_el/(hbar*Vf),1.5*e_el/(hbar*Vf),100);

% y1 = hbar*Vf*x/e_el+Vi;
% y2 = -hbar*Vf*x/e_el+Vi;
 y1 = hbar*Vf*x/e_el;
 y2 = -hbar*Vf*x/e_el;

%x1 = 2*x(end) + linspace(-1.5*max(abs(ktop+k_i),abs(kbott+k_i)),1.5*max(abs(ktop+k_i),abs(kbott+k_i)),100);
x1 = 2*x(end) + x;
% 
% y3 = hbar*Vf*(x1-2*x(end))/e_el;%+Vi;
% y4 = -hbar*Vf*(x1-2*x(end))/e_el;%+Vi;

y3 = hbar*Vf*(x1-2*x(end))/e_el-Vi; %like adding V_kp
y4 = -hbar*Vf*(x1-2*x(end))/e_el-Vi;
if flag ==0
figure;
end

plot(x,y1,'b');
hold on
plot(x,y2,'b');
plot(x1,y3,'r');
plot(x1,y4,'r');
%yline(Ef_bott+Vi,'b')
%yline(Ef_bott,'b')
yplot_value(x,Ef_bott,'b')

%yline(Ef_top,'r')
%yline(Ef_top-Vi,'r')
yplot_value(x1,Ef_top-Vi,'r')
%yline(Vi/2,'--g')
yline(-Vi/2,'--g')
title(['Bandstructure at V_b =' ,num2str(V,'%.2f'),' V and V_g =',num2str(Vg,'%.2f'),' V'])
hold off

end