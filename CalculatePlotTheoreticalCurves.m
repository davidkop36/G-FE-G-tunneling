function [V_bias_zero_Q_bottom,V_gate_zero_Q_top,V_gate_zero_Vi]= CalculatePlotTheoreticalCurves(V_all,Vg_all,V_kp_0,C0,C1,beta,er,Qbottom_zero_locus,Qtop_zero_locus,Vi_zero_locus,theta)


V_bias_zero_Q_bottom =-V_kp_0+C1/(C0*er)*Vg_all+beta*sqrt(C1)*sqrt(abs(Vg_all)).*sign(Vg_all);
figure;
plot(V_all(Qbottom_zero_locus),Vg_all,'-','LineWidth',2);
hold on
plot(V_bias_zero_Q_bottom,Vg_all,'--','LineWidth',2);
xline(V_all(1));
xline(V_all(end));
legend('Numerical curve','Theoretical curve','Lower cutoff in numerical','Upper cutoff in numerical','Location','bestoutside')
xlabel('V_{bias}');
ylabel('V_{gate}');
title('Locus of points where Q_{bottom} is zero')

V_gate_zero_Q_top = (V_all+V_kp_0).* (-abs( V_all+V_kp_0) /(beta^2*C1)-1  );
figure;
plot(V_all(Qtop_zero_locus),Vg_all,'-','LineWidth',2);
hold on
plot(V_all,V_gate_zero_Q_top,'--','LineWidth',2);
yline(Vg_all(1));
yline(Vg_all(end));
legend('Numerical curve','Theoretical curve','Lower cutoff in numerical','Upper cutoff in numerical','Location','bestoutside')
xlabel('V_{bias}');
ylabel('V_{gate}');
title('Locus of points where Q_{top} is zero')

if theta==0
V_gate_zero_Vi= -V_all + beta*sqrt(abs(er * C0*V_kp_0)).*sign(V_kp_0) + er*(C0/C1)*V_kp_0 +(1/C1)* ((-V_all/(beta)  + sqrt(abs(er *C0*V_kp_0)).*sign(V_kp_0) ).^2).*sign((-V_all/(beta)  + sqrt(abs(er *C0*V_kp_0)).*sign(V_kp_0) ) ) ;
figure;
plot(V_all(Vi_zero_locus),Vg_all,'-','LineWidth',2);
hold on
plot(V_all,V_gate_zero_Vi,'--','LineWidth',2);
yline(Vg_all(1));
yline(Vg_all(end));
legend('Numerical curve','Theoretical curve','Lower cutoff in numerical','Upper cutoff in numerical','Location','bestoutside')
xlabel('V_{bias}');
ylabel('V_{gate}');
title('Locus of points where V_i is zero')
else
    V_gate_zero_Vi=0;
end

end