clear all

% Input parameters of the model
er = 3; %FE material dielectric constant
er_gate = 3.8; %dielectric constant of the gate material
d = 5e-9; % Junction thickness. Default to 5nm . Try 1nm, per Moshe's recommendation 
d_gate = 90e-9; %gate thickness, 15e-9 standard 
theta  = 2 *pi/180; %twist angle in between two graphene sheets
V_kp_0 =+0.1 % Bare ferroelctricity voltage

% Bias_gate range of interest
V_b_max = 2;
V_b_min = -2;
V_g_max = 30;
V_g_min = -30;

course_grid = [101, 103]; %used for electrostatics calculation
fine_grid = [1001, 1003]; %used for I-V current calculation


%Physical constants
e0 = 8.854e-12;
e_el = 1.602e-19;
Vf = 1e6; % Fermi velocity of graphene electrones
hbar = 1.055e-34;


%% Execution part

% Derived parameters of the model
beta = sqrt(pi/e_el^3)*hbar*Vf; 
C0 = e0/d;
C1 = e0*er_gate/d_gate;

V_all = linspace(V_b_min,V_b_max,course_grid(1)); 
Vg_all = linspace(V_g_min,V_g_max,course_grid(2)); 


%Custom functions 
symsq = @(x) sqrt(abs(x))*sign(x); 
syms Q Q1


for ii=1:length(V_all)
    V = V_all(ii);
    for jj=1:length(Vg_all)
        % Solve for positive polarization
        Vg = Vg_all(jj);
        eq1 = (Q/(C0*er) - V - V_kp_0 + beta * (symsq(Q) - symsq(Q1-Q))) == 0;
        eq2  = (Q1/C1 - Vg + beta * (symsq(Q1-Q))) == 0;
        [Q_eq,Q1_eq] = vpasolve([eq1,eq2],[Q,Q1]);
        Q_all(ii,jj) = Q_eq;
        Q1_all(ii,jj) = Q1_eq;
        
        % Solve for negative polarization
        eq1 = (Q/(C0*er) - V + V_kp_0 + beta * (symsq(Q) - symsq(Q1-Q))) == 0;
        eq2  = (Q1/C1 - Vg + beta * (symsq(Q1-Q))) == 0;
        [Q_eq,Q1_eq] = vpasolve([eq1,eq2],[Q,Q1]);
        Q_all2(ii,jj) = Q_eq;
        Q1_all2(ii,jj) = Q1_eq;
        
        % Solve for the absence of the polarization
        eq1 = (Q/(C0*er) - V  + beta * (symsq(Q) - symsq(Q1-Q))) == 0;
        eq2  = (Q1/C1 - Vg + beta * (symsq(Q1-Q))) == 0;
        [Q_eq3,Q1_eq3] = vpasolve([eq1,eq2],[Q,Q1]);
        Q_all3(ii,jj) = Q_eq3;
        Q1_all3(ii,jj) = Q1_eq3;   
        
    end
end

Q_all = eval(Q_all);
Q1_all = eval(Q1_all);

Q_all2 = eval(Q_all2);
Q1_all2 = eval(Q1_all2);

Q_all3 = eval(Q_all3);
Q1_all3 = eval(Q1_all3);

%% Interpolate to higher resolution
Vg_fine = linspace(Vg_all(1),  Vg_all(end),fine_grid(2));
V_fine = linspace(V_all(1),  V_all(end),fine_grid(1));
[xq,yq] = ndgrid(V_fine,Vg_fine);
[x0,y0] = ndgrid(V_all,Vg_all);

F = griddedInterpolant(x0,y0,Q_all);
Q_all = F(xq,yq);
F = griddedInterpolant(x0,y0,Q1_all);
Q1_all = F(xq,yq);

F2 = griddedInterpolant(x0,y0,Q_all2);
Q_all2 = F2(xq,yq);
F2 = griddedInterpolant(x0,y0,Q1_all2);
Q1_all2 = F2(xq,yq);


F3 = griddedInterpolant(x0,y0,Q_all3);
Q_all3 = F3(xq,yq);
F3 = griddedInterpolant(x0,y0,Q1_all3);
Q1_all3 = F3(xq,yq);


V_all = V_fine;
Vg_all = Vg_fine;

% Calculate the energetical mismatch between the Dirac cones 
Vkp = - Q_all/(C0*er) + V_kp_0;
Vkp2 = - Q_all2/(C0*er) - V_kp_0;
Vkp3 = - Q_all3/(C0*er);

%Visulize the difference in KP voltage
figure;
f_thing = (Vkp-Vkp2)/(2*V_kp_0);
contourf(V_all,Vg_all,f_thing',12)
ylabel('Gate voltage [V]','FontSize',15);
xlabel('Bias voltage [V]','FontSize',15);
set(gca,'YDir','normal')
ax = gca;
ax.FontSize = 20;
colorbar;

%% Preview charge plots
figure;
imagesc(V_all,Vg_all,(Q1_all-Q_all));
xlabel('V_b');
ylabel('V_g');
colorbar;
title('Charge on the bottom graphene');
set(gca,'YDir','normal');
figure;
imagesc(V_all,Vg_all,(Q_all)');
xlabel('V_b');
ylabel('V_g');
title('Charge on the top graphene')
colorbar;
set(gca,'YDir','normal');


%%
ind=600;
figure;
diffs = (Vkp-Vkp2)./(2*V_kp_0);
diffs_zero = diffs(ind,:);
plot(Vg_all,diffs(ind,:),'LineWidth',2)
xlabel('Gate voltage [V]','FontSize',20)
ylabel(' $\frac{V_{kp}^+ - V_{kp}^-}{2*V_{kp}^0}$','FontSize',20,'interpreter','latex')
title(['For Vb = ',num2str(V_all(ind))  ],'FontSize',20)
aa = 0.5*sqrt(beta*er*C0/sqrt(C1));
g_plus = V_all(ind) + V_kp_0 +beta*sqrt(abs(Vg_all)*C1).*sign(Vg_all);
g_minus= V_all(ind) - V_kp_0 +beta*sqrt(abs(Vg_all)*C1).*sign(Vg_all);

asymp = aa*(abs(Vg_all)).^(-1/4);
asymp2 = 2*beta./(sqrt(beta^2+4*abs(g_plus)/(er*C0))  + sqrt(beta^2+4*abs(g_minus)/(er*C0))  );

hold on
plot(Vg_all(Vg_all>0.5),asymp(Vg_all > 0.5),'LineWidth',2)
legend('Numerical curve','Asymptotic expansion','Location','bestoutside','FontSize',15)
%plot(Vg_all(Vg_all>0),asymp2(Vg_all > 0))
ax = gca;
ax.FontSize = 20;

hold off
%% Cut along fixed gate

indg=700;
figure;
diffs = (Vkp-Vkp2)./(2*V_kp_0);
diffs_zero = diffs(:,indg);
plot(V_all,diffs(:,indg),'LineWidth',2);
xlabel('Bias voltage [V]','FontSize',20)
ylabel(' $\frac{V_{kp}^+ - V_{kp}^-}{2*V_{kp}^0}$','FontSize',20,'interpreter','latex')
title(['For Vg = ',num2str(Vg_all(indg),'%.2f'),' V'  ],'FontSize',20)
aa = 0.5*sqrt(beta*er*C0/sqrt(C1));
g_plus = V_all+ V_kp_0 +beta*sqrt(abs(Vg_all(indg))*C1).*sign(Vg_all(indg));
g_minus= V_all - V_kp_0 +beta*sqrt(abs(Vg_all(indg))*C1).*sign(Vg_all(indg));

%asymp = aa*(abs(Vg_all)).^(-1/4);
%asymp2 = 2*beta./(sqrt(beta^2+4*abs(g_plus)/(er*C0))  + sqrt(beta^2+4*abs(g_minus)/(er*C0))  );

legend('Numerical curve','Asymptotic expansion','Location','bestoutside','FontSize',15)
%plot(Vg_all(Vg_all>0),asymp2(Vg_all > 0))
ax = gca;
ax.FontSize = 20;

hold off


%%
indg =600;
Q0 = (V_all+V_kp_0)*er*C0;
Q0_prim = Vg_all(indg)*C1;
% a1 = 1/(er*C0);
% a2 = 1/C1;
% b1= beta./(2*sqrt(abs(Q0)));
% b2= beta./(2*sqrt(abs(Q0_prim-Q0)));
% c1 = -beta*sqrt(abs(Q0)).*sign(Q0);
% c2 =  beta*(2*sqrt(abs(Q0_prim-Q0))).*sign(Q0_prim-Q0);
% 
% delta_Q  = (a2.*c1+b2.*c1+a2.*c2)./(a1.*a2 + a1.*b2 + a2.*b1+b1.*b2+a2.*b2);
% delta_Q_prim  = (b2.*c1-a1.*c2-b1.*c2)./(a1.*a2 + a1.*b2 + a2.*b1+b1.*b2+a2.*b2);
% 
% Q02 = (V_all - V_kp_0)*er*C0;
% a1 = 1/(er*C0);
% a2 = 1/C1;
% b1= beta./(2*sqrt(abs(Q02)));
% b2= beta./(2*sqrt(abs(Q0_prim-Q02)));
% c1 = -beta*sqrt(abs(Q02)).*sign(Q02);
% c2 =  beta*(2*sqrt(abs(Q0_prim-Q02))).*sign(Q0_prim-Q02);
% 
% delta_Q2  = (a2.*c1+b2.*c1+a2.*c2)./(a1.*a2 + a1.*b2 + a2.*b1+b1.*b2+a2.*b2);
% delta_Q_prim2  = (b2.*c1-a1.*c2-b1.*c2)./(a1.*a2 + a1.*b2 + a2.*b1+b1.*b2+a2.*b2);
% 
% 
% f_fit = 1 - (Q0+delta_Q-Q02-delta_Q2)./(2*V_kp_0*er*C0);
%f_fit2 =  - (delta_Q-delta_Q2)./(2*V_kp_0*er*C0);
f_fit2 =  beta*sqrt(er*C0./V_all);

figure;
diffs = (Vkp-Vkp2)./(2*V_kp_0);
diffs_zero = diffs(:,indg);
plot(V_all,diffs(:,indg),'LineWidth',2);
xlabel('V_b [V]','FontSize',20)
ylabel(' f','FontSize',20)
title(['For V_g = ',num2str(Vg_all(indg),'%.2f'),' V'  ],'FontSize',20)
hold on
plot(V_all((V_all)>V_kp_0*4),f_fit2((V_all)>V_kp_0*4),'r--','LineWidth',2)

f_fit_peak = 2./(sqrt(1+4*abs(V_all+V_kp_0)./(C0*er*beta^2)) +sqrt(1+4*abs(V_all-V_kp_0)./(C0*er*beta^2))  );

%plot(V_all(-V_all*sign(Vg_all(indg)) > 0),f_fit_peak(-V_all*sign(Vg_all(indg)) > 0),'g.-','LineWidth',2)


%%
% figure
% plot(Q1_all(ind,:))
% hold on
% plot(Q1_all2(ind,:))
% 
% hold off

% In terms of charge
% figure
% plot(Q1_all(ind,:),Vkp(ind,:)/V_kp_0)
% hold on
% plot(Q1_all2(ind,:),Vkp2(ind,:)/V_kp_0)
% hold off
Q_bottom= Q1_all-Q_all;
Q_bottom_2 = Q1_all2-Q_all2;
Q_bottom_3 = Q1_all3-Q_all3;

figure;
Vkp2_intp = interp1(Q_bottom_2(ind,:),Vkp2(ind,:),Q_bottom(ind,:));
plot(Q_bottom(ind,:),(Vkp(ind,:)-Vkp2_intp)./(2*V_kp_0))  ;
xlabel('Q_{bottom}')
ylabel('(Vkp-Vkp2)./(2*V_{kp}^0)')
title(['For Vb = ',num2str(V_all(ind))  ])
hold on




%% Find Fermi energies relative to Dirac point
% up
Ef_top = -sqrt(4*pi*abs(Q_all)/e_el).*sign(Q_all)*hbar*Vf/(2*e_el);
Ef_bott = -sqrt(4*pi*abs(Q_bottom)/e_el).*sign(Q_bottom)*hbar*Vf/(2*e_el);
V_all_rep = repmat(V_all',1,length(Vg_all));
Vi = V_all_rep + Ef_top - Ef_bott;
%Vi = V_all_rep - Ef_top + Ef_bott;

qc = (12e-9)^(-1);
gamma = e_el*(Ef_bott-Ef_top-V_all_rep)/(hbar*Vf*qc);
kt1 = Ef_top*e_el/(hbar*Vf*qc);
%kt2 = (Ef_top+V_all_rep)*e_el/(hbar*Vf*qc);
kt2 = (Ef_top+V_all_rep)*e_el/(hbar*Vf*qc);
k_vi = e_el*Vi/(hbar*Vf*qc);
%k_vi = -e_el*Vi/(hbar*Vf*qc);

%down
Ef_top_2 =- sqrt(4*pi*abs(Q_all2)/e_el).*sign(Q_all2)*hbar*Vf/(2*e_el);
Ef_bott_2 = -sqrt(4*pi*abs(Q_bottom_2)/e_el).*sign(Q_bottom_2)*hbar*Vf/(2*e_el);
Ef_top_3 =- sqrt(4*pi*abs(Q_all3)/e_el).*sign(Q_all3)*hbar*Vf/(2*e_el);
Ef_bott_3 =- sqrt(4*pi*abs(Q_bottom_3)/e_el).*sign(Q_bottom_3)*hbar*Vf/(2*e_el);
%V_all_rep = repmat(V_all',1,length(Vg_all));
Vi_2 = V_all_rep + Ef_top_2 - Ef_bott_2;
Vi_3 = V_all_rep + Ef_top_3 - Ef_bott_3;

%Vi_2 = V_all_rep - Ef_top_2 + Ef_bott_2;
%Vi_3 = V_all_rep - Ef_top_3 + Ef_bott_3;

gamma_2 = e_el*(Ef_bott_2-Ef_top_2-V_all_rep)/(hbar*Vf*qc);
kt1_2 = Ef_top_2*e_el/(hbar*Vf*qc);
%kt2_2 = (Ef_top_2+V_all_rep)*e_el/(hbar*Vf*qc);
kt2_2 = (Ef_top_2+V_all_rep)*e_el/(hbar*Vf*qc);

k_vi_2 = e_el*Vi_2/(hbar*Vf*qc);
%k_vi_2 = -e_el*Vi_2/(hbar*Vf*qc);

% for ii=1:size(kt1,1)
%     for jj=1:size(kt1,2)
%         k = linspace(kt1(ii,jj),kt2(ii,jj),1000);
% 
%         I1(ii,jj) = trapz(k, abs(k-k_vi(ii,jj)).*abs(k).*(1+k.^2 + (k - k_vi(ii,jj) ).^2  ) ./ ( (1 + k_vi(ii,jj)^2)^(3/2) .* (1 + (2*k-k_vi(ii,jj)).^2 ).^(3/2)    ));
% 
% 
%     end
% end

%I1 = Calculate_current_graphene(kt1,kt2,k_vi);
%I2 = Calculate_current_graphene(kt1_2,kt2_2,k_vi_2);
reduced_twist_vector = (8*pi/(3*1.42e-10))*sin(theta/2)/qc;
I1 = Calculate_current_graphene(kt1,kt2,k_vi,reduced_twist_vector);
I2 = Calculate_current_graphene(kt1_2,kt2_2,k_vi_2,reduced_twist_vector);

gde = 502;
locs_all_local_max = find( islocalmax(abs(I1(:,gde))) == 1);
max_location_peak = locs_all_local_max(find(max(abs(I1(locs_all_local_max,gde))) == abs(I1(locs_all_local_max,gde))  ));
figure;
plot(V_all,I1(:,gde),'LineWidth',2)
%
peak = (1./(1 + k_vi(:,gde).^2)).^(3/2);
peak = (peak./max(peak)).*I1(max_location_peak,gde);
hold on
plot(V_all,peak,'LineWidth',2 )
hold off
xlabel('Bias [V]','FontSize',20)
ylabel('Current [a.u.]','FontSize',20)
title(['Current vs. bias for Vg =  ',num2str(Vg_all(gde))] ,'FontSize',20)
legend('Numerical I-V curve','Prefactor contribution to the peak','Location','bestoutside','FontSize',15)


figure;
surf(V_all,Vg_all,real(I1)')
xlabel('Bias [V]')
ylabel('Gate voltage [V]')
zlabel('Current [a.u.]')

figure;
contour(V_all,Vg_all,real(I1)',20)
xlabel('Bias [V]')
ylabel('Gate voltage [V]')

% figure;
% contour(V_all,Vg_all,real(I1-I2)',20)
% xlabel('Bias [V]')
% ylabel('Gate voltage [V]')
figure;
imagesc(V_all,Vg_all,abs(I1)')
xlabel('Bias [V]')
ylabel('Gate voltage [V]')





%%
gde1 =550;

figure;
plot(V_all,I1(:,gde1),'LineWidth',2)
hold on
plot(V_all,I2(:,gde1),'LineWidth',2)
xlabel('Bias [V]','FontSize',20)
ylabel('Current [a.u.]','FontSize',20)
title(['Current vs. bias for Vg =  ',num2str(Vg_all(gde1),'%.2f'),'V'],'FontSize',20 )
legend('V_{KP}^{(0)} > 0 ','V_{KP}^{(0)} < 0 ' ,'Location','bestoutside','FontSize',15)
ax = gca;
ax.FontSize = 15
%%
figure;
imagesc(V_all,Vg_all,abs(Vi'))
figure;
surf(V_all,Vg_all,abs(Q_all'))

figure;
imagesc(V_all,Vg_all,abs(I1'))
hold on
%Vi_zero_locus = ExtractMinimumInDirection(abs(Vi),1);
 Vi_zero_locus = ExtractMinimumInDirection(abs(Vi-reduced_twist_vector*qc*hbar*Vf/e_el),1);
  Vi_zero_locus_other = ExtractMinimumInDirection(abs(Vi+reduced_twist_vector*qc*hbar*Vf/e_el),1);

Qtop_zero_locus = ExtractMinimumInDirection(abs(Q_all),1);
Qbottom_zero_locus = ExtractMinimumInDirection(abs(Q_bottom),1);

plot(V_all(Vi_zero_locus),Vg_all,'r')
plot(V_all(Qtop_zero_locus),Vg_all,'g')
plot(V_all(Qbottom_zero_locus),Vg_all,'y')

hold off
%%




Vi_zero_locus_2 = ExtractMinimumInDirection(abs(Vi_2-reduced_twist_vector*qc*hbar*Vf/e_el),1);
Vi_zero_locus_other_2= ExtractMinimumInDirection(abs(Vi_2+reduced_twist_vector*qc*hbar*Vf/e_el),1);
Qtop_zero_locus_2 = ExtractMinimumInDirection(abs(Q_all2),1);
Qbottom_zero_locus_2 = ExtractMinimumInDirection(abs(Q_bottom_2),1);
 Vi_zero_locus_0 = ExtractMinimumInDirection(abs(Vi_3-reduced_twist_vector*qc*hbar*Vf/e_el),1);
 Qtop_zero_locus_0 = ExtractMinimumInDirection(abs(Q_all3),1);
 Qbottom_zero_locus_0 = ExtractMinimumInDirection(abs(Q_bottom_3),1);
 Qprim_zero_locus_0 = ExtractMinimumInDirection(abs(Q1_all),1);



Qcompens_zero_locus_0 = ExtractMinimumInDirection(abs(Q_all + Q_all2),1);


[V_bias_zero_Q_bottom,V_gate_zero_Q_top,V_gate_zero_Vi]= CalculatePlotTheoreticalCurves(V_all,Vg_all,V_kp_0,C0,C1,beta,er,Qbottom_zero_locus,Qtop_zero_locus,Vi_zero_locus,theta);
[V_bias_zero_Q_bottom2,V_gate_zero_Q_top2,V_gate_zero_Vi2]= CalculatePlotTheoreticalCurves(V_all,Vg_all,-V_kp_0,C0,C1,beta,er,Qbottom_zero_locus_2,Qtop_zero_locus_2,Vi_zero_locus_2,theta);
%[V_bias_zero_Q_bottom0,V_gate_zero_Q_top0,V_gate_zero_Vi0]= CalculatePlotTheoreticalCurves(V_all,Vg_all,0,C0,C1,beta,er,Qbottom_zero_locus_0,Qtop_zero_locus_0,Vi_zero_locus_0);


% Theoretical curves over
%%

plot_type_option = 1;

figure;
if plot_type_option==2
%imagesc(V_all,Vg_all,(I2')+(I1'))
uslov = abs(V_all) < 2;
imagesc(V_all(uslov),Vg_all,(I2(uslov,:)')+(I1(uslov,:)'))

else
imagesc(V_all,Vg_all,(I1'))

end
set(gca,'YDir','normal')


hold on
lw_def = 2;

% 
% 
% figure;
% imagesc(real(I_twist)')
% xlabel('Bias [V]')
% ylabel('Gate voltage [V]')
% set(gca,'YDir','normal')


% plot(V_all(Vi_zero_locus),Vg_all,'r--','LineWidth',lw_def)
% plot(V_all(Qtop_zero_locus),Vg_all,'g--','LineWidth',lw_def)
% plot(V_all(Qbottom_zero_locus),Vg_all,'y--','LineWidth',lw_def)
% plot(V_all(Vi_zero_locus_2),Vg_all,'r-.','LineWidth',lw_def)
% plot(V_all(Qtop_zero_locus_2),Vg_all,'g-.','LineWidth',lw_def)
% plot(V_all(Qbottom_zero_locus_2),Vg_all,'y-.','LineWidth',lw_def)
% %plot(V_all(Vi_zero_locus_0),Vg_all,'r','LineWidth',lw_def)
% plot(V_all(Qtop_zero_locus_0),Vg_all,'g','LineWidth',lw_def)
% plot(V_all(Qbottom_zero_locus_0),Vg_all,'y','LineWidth',lw_def)  


%theoretical curves
plot(V_all,V_gate_zero_Q_top,'g--','LineWidth',lw_def)
plot(V_bias_zero_Q_bottom,Vg_all,'y--','LineWidth',lw_def)
%plot(V_all,V_gate_zero_Vi,'r--','LineWidth',lw_def)
if theta ==0
plot(V_all,Vi_zero_locus,'r--','LineWidth',lw_def)
else
plot(V_all(Vi_zero_locus),Vg_all,'r--','LineWidth',lw_def)
plot(V_all(Vi_zero_locus_other),Vg_all,'r--','LineWidth',lw_def)

end    

plot(V_all,V_gate_zero_Q_top2,'g.','LineWidth',lw_def)
plot(V_bias_zero_Q_bottom2,Vg_all,'y.','LineWidth',lw_def)
%plot(V_all,V_gate_zero_Vi2,'r.','LineWidth',lw_def)
if theta==0
plot(V_all,Vi_zero_locus_2,'r--','LineWidth',lw_def)
else
plot(V_all(Vi_zero_locus_2),Vg_all,'r.','LineWidth',lw_def)
plot(V_all(Vi_zero_locus_other_2),Vg_all,'r.','LineWidth',lw_def)


end 

%plot(V_all(Qcompens_zero_locus_0),Vg_all,'m','LineWidth',lw_def)  
%plot(V_all,V_gate_zero_Vi0,'b--','LineWidth',lw_def)  
%plot(V_all(Vi_zero_locus_0),Vg_all,'b','LineWidth',lw_def)




xlabel('Bias [V]','FontSize', 20)
ylabel('Gate [V]','FontSize', 20)
%%legend('V_i = 0, P>0',' Q_{top} = 0, P>0','Q_{bottom} = 0, P>0','V_i = 0, P<0', 'Q_{top} = 0, P<0','Q_{bottom} = 0, P<0','V_i = 0, P=0',' Q_{top}=0, P=0',' Q_{bottom}=0, P=0','Location','bestoutside'	)
%legend('V_i = 0, P>0',' Q_{top} = 0, P>0','Q_{bottom} = 0, P>0','V_i = 0, P<0', 'Q_{top} = 0, P<0','Q_{bottom} = 0, P<0',' Q_{top}=0, P=0',' Q_{bottom}=0, P=0','Location','bestoutside','FontSize', 15	)
if theta ==0
    legend('Q_{top}^{+}=0',' Q_{bottom}^{+}=0',' V_{KP}^{+}=0','Q_{top}^{-}=0',' Q_{bottom}^{-}=0',' V_{KP}^{-}=0','Location','bestoutside','FontSize', 15	)
else
    legend('Q_{top}^{+}=0',' Q_{bottom}^{+}=0',' V_{KP}^{+}=\Delta',' V_{KP}^{-}=\Delta','Q_{top}^{-}=0',' Q_{bottom}^{-}=0',' V_{KP}^{+}=-\Delta',' V_{KP}^{-}=-\Delta','Location','bestoutside','FontSize', 15	)

end

ax = gca;
ax.FontSize = 15;
title('I_{+}','FontSize', 20)
colorbar;
%% Plot at fixed gate
closestIndex2 = 836; %Fix at 20V gate
midindex = 661; %Fix at 9.52V gate

%midindex = 502;
yline(Vg_all(closestIndex2));
yline(Vg_all(midindex));


if theta == 0
[x0,y0] = ginput(6);
scatter(x0,repmat(Vg_all(closestIndex2),1,6))
scatter(x0,repmat(Vg_all(midindex),1,6))
else
[x0,y0] = ginput(4);
scatter(x0,repmat(Vg_all(closestIndex2),1,4))
scatter(x0,repmat(Vg_all(midindex),1,4)) %comment out later

end
hold off
figure;
for ii=1:4
[minValue,closestIndex] = min(abs(V_all-x0(ii)));

if theta == 0
graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi(closestIndex,closestIndex2),Q_all(closestIndex,closestIndex2),Q_bottom(closestIndex,closestIndex2));
graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(midindex),Vi(closestIndex,midindex),Q_all(closestIndex,midindex),Q_bottom(closestIndex,midindex));
% 2nd line is the zero bias plot
else
graphene_plot_band_structure_qc_twisted_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi(closestIndex,closestIndex2),Q_all(closestIndex,closestIndex2),Q_bottom(closestIndex,closestIndex2),theta);
%comment out the 2nd line later
end
end


%%
figure;

plot(V_all,I1(:,closestIndex2)','-b','LineWidth',lw_def)
hold on
plot(V_all,I2(:,closestIndex2)','--r','LineWidth',lw_def)
xlabel('V_b [V]','FontSize', 20)
ylabel('I [a.u.]','FontSize', 20)
legend('I^{(+)}','I^{(-)}','Location','bestoutside','FontSize', 15	)
title(['I-V curve for ', num2str(Vg_all(closestIndex2),'%.2f'), ' V'],'FontSize', 20)

figure;

plot(V_all,I1(:,midindex)','-b','LineWidth',lw_def)
hold on
plot(V_all,I2(:,midindex)','--r','LineWidth',lw_def)
xlabel('V_b [V]','FontSize', 20)
ylabel('I [a.u.]','FontSize', 20)
legend('I^{(+)}','I^{(-)}','Location','bestoutside','FontSize', 15	)
title(['I-V curve for ', num2str(Vg_all(midindex),'%.2f'), ' V'],'FontSize', 20)

%%
% Debug
closestIndex2 = 836; %Fix at 20V gate
yline(Vg_all(closestIndex2));
[x0,y0] = ginput(5);
scatter(x0,y0)
hold off
figure;
for ii=1:5
[minValue,closestIndex] = min(abs(V_all-x0(ii)));
[minValue2,closestIndex2] = min(abs(Vg_all-y0(ii)));
graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi(closestIndex,closestIndex2),Q_all(closestIndex,closestIndex2),Q_bottom(closestIndex,closestIndex2));

end

% closestIndex = 860;
% 
% graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi(closestIndex,closestIndex2),Q_all(closestIndex,closestIndex2),Q_bottom(closestIndex,closestIndex2));
% graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi_2(closestIndex,closestIndex2),Q_all2(closestIndex,closestIndex2),Q_bottom_2(closestIndex,closestIndex2));
% graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi_3(closestIndex,closestIndex2),Q_all3(closestIndex,closestIndex2),Q_bottom_3(closestIndex,closestIndex2));



%graphene_plot_band_structure_qc_debug(V_all(closestIndex),Vg_all(closestIndex2),0,0,0);


%% Distance between the peaks
figure
plot(Vg_all,V_all(Vi_zero_locus) -V_all(Vi_zero_locus_2),'LineWidth',lw_def)
diff_theor = 2*beta*(er*C0*V_kp_0).^(1/2) + 2*beta*er*C0*V_kp_0./sqrt(beta^2*C1^2-4*Vg_all*C1);
hold on
plot(Vg_all(Vg_all<0),diff_theor(Vg_all<0),'LineWidth',lw_def   )
xlabel('Gate voltage [V]','Fontsize',20)
ylabel('Difference in peak positions [V]','Fontsize',20)
legend('Numerical result','Asymptotic dependence derived from EOM-s','Fontsize',15)
ax = gca;
ax.FontSize = 15;


%% Debug



%


%% Saving

I_V.V_all= V_all;
I_V.Vg_all= Vg_all;
I_V.I1=I1;
I_V.I2=I2;

I_V.Q_all=Q_all;
I_V.Q_bottom=Q_bottom;
I_V.Q_all_2=Q_all2;
I_V.Q_bottom_2=Q_bottom_2;
I_V.Q_all_3=Q_all3;
I_V.Q_bottom_3=Q_bottom_3;
I_V.Vi=Vi;
I_V.Vi2=Vi_2;
I_V.Vi0=Vi_3;
I_V.qc=qc;
I_V.Vkp=Vkp;
I_V.Vkp2=Vkp2;
I_V.Vkp3=Vkp3;


save('I_V_16','I_V')
%[V_all,Vg_all,I1,Q_all,Q_bottom,Vi,qc,Vkp] =   ExtractSavedParamters(I_V);
%%
% 
% xx  = Q_bottom(ind,:);
% yy = (Vkp(ind,:)-Vkp2_intp)./(2*V_kp_0);
% 
% 
% %%
% kbT = 0.025;
% q_max = 0.1*  4*pi/(3*0.142);
% z_faktor = 500;
% kb = linspace(0,q_max,1000);
% kt = linspace(0,q_max,1000);
% %angles = linspace(0,2*pi,1000);
% %Vf = 1e6;
% hbar2 =6.582e-16 * 1e9;
% %kbt = 0.025;
% % zero gate
% 
% 
% for ii=1:length(V_all)
%     for jj=1:length(Vg_all)
% mu_b = Ef_bott(ii,jj); % +P
% mu_t = Ef_top(ii,jj); % +P
% % mu_b_1 = Vg - deltaE_all_1(ind,ii); %-P
% % mu_t_1 = Vg + deltaE_all_1(ind,ii); % -P
% 
% fb_qb = 1./(1+exp((kb*Vf*hbar2-mu_b)/kbT));
% ft_qb = 1./(1+exp((kt*Vf*hbar2-mu_t)/kbT));
% ft_minus_qb = 1./(1+exp((-kb*Vf*hbar2-mu_b)/kbT));
% ft_minus_qt = 1./(1+exp((-kt*Vf*hbar2-mu_t)/kbT));
% 
% % fb_qb_1 = 1./(1+exp((kb*Vf*hbar2-mu_b_1)/kbT)); %-P
% % ft_qb_1 = 1./(1+exp((kt*Vf*hbar2-mu_t_1)/kbT));
% % ft_minus_qb_1 = 1./(1+exp((-kb*Vf*hbar2-mu_b_1)/kbT));
% % ft_minus_qt_1 = 1./(1+exp((-kt*Vf*hbar2-mu_t_1)/kbT));
% 
% j = zeros(length(kb),length(kt));
% j2 = j;
% j3 = j;
% 
% j = -abs(kb'*kt).*(qc^2+(kb').^2+ kt.^2).*(1./(qc^2+(kb'+kt).^2).^(3/2) .* 1./(qc^2+(kb'-kt).^2)^(3/2) )  .* (abs((kb'-kt)*hbar2*Vf+Vi(ind,ii)) < q_max/z_faktor ) .* (fb_qb' -ft_qb) ;% + (1/qc^2-1./(qc^2+4*kb.^2)).*(fb_qb-ft_minus_qb);
% j2 = -abs(kb'*kt).*(qc^2+(kb').^2+ kt.^2).*(1./(qc^2+(kb'+kt).^2).^(3/2) .* 1./(qc^2+(kb'-kt).^2)^(3/2))  .* (abs((kb'-kt)*hbar2*Vf-Vi(ind,ii)) < q_max/z_faktor ) .* (ft_minus_qb' -ft_minus_qt);
% %hole hole
% %j1_1 = -(1./(qc^2+(kb'+kt).^2) - 1./(qc^2+(kb'-kt).^2))  .* (abs((kb'-kt)*hbar2*Vf+Vi_1(ind,ii)) < q_max/z_faktor ) .* (fb_qb_1' -ft_qb_1) ;% + (1/qc^2-1./(qc^2+4*kb.^2)).*(fb_qb-ft_minus_qb);
% %j2_1 = -(1./(qc^2+(kb'+kt).^2) - 1./(qc^2+(kb'-kt).^2))  .* (abs((kb'-kt)*hbar2*Vf-Vi_1(ind,ii)) < q_max/z_faktor ) .* (ft_minus_qb_1' -ft_minus_qt_1);
% 
% if Vi(ii,jj)>=0
% %j3  =- (1./(qc^2+(kb'+kt).^2) - 1./(qc^2+(kb'-kt).^2)).* (abs((kb'+kt)*hbar2*Vf-Vi(ind,ii)) < q_max/1000 ) .* (fb_qb' -ft_minus_qt);  % electron-hole
% j3  = -abs(kb'*kt).*(qc^2+(kb').^2+ kt.^2).*(1./(qc^2+(kb'+kt).^2).^(3/2) .* 1./(qc^2+(kb'-kt).^2).^(3/2)).* (abs((kb'+kt)*hbar2*Vf-Vi(ind,ii)) < q_max/z_faktor ) .* (ft_minus_qb' -ft_qb);  % electron-hole    
% 
% else
% %j3  = -(1./(qc^2+(kb'+kt).^2) - 1./(qc^2+(kb'-kt).^2)).* (abs((kb'+kt)*hbar2*Vf-Vi(ind,ii)) < q_max/1000 ) .* (ft_minus_qb' -ft_qb);  % electron-hole    
% j3  = -abs(kb'*kt).*(qc^2+(kb').^2+ kt.^2).*(1./(qc^2+(kb'+kt).^2).^(3/2) .* 1./(qc^2+(kb'-kt).^2).^(3/2)).* (abs((kb'+kt)*hbar2*Vf+Vi(ind,ii)) < q_max/z_faktor ) .* (fb_qb' -ft_minus_qt);  % electron-hole
% end
% 
% % if Vi_1(ind,ii)>=0
% % j3_1  = -(1./(qc^2+(kb'+kt).^2) - 1./(qc^2+(kb'-kt).^2)).* (abs((kb'+kt)*hbar2*Vf-Vi_1(ind,ii)) < q_max/z_faktor ) .* (ft_minus_qb_1' -ft_qb_1);  % electron-hole    
% % else
% %  j3_1  =- (1./(qc^2+(kb'+kt).^2) - 1./(qc^2+(kb'-kt).^2)).* (abs((kb'+kt)*hbar2*Vf+Vi_1(ind,ii)) < q_max/z_faktor ) .* (fb_qb_1' -ft_minus_qt_1);  % electron-hole
% % end   
% 
% 
% 
% 
% j_all_1(ii,jj) = sum(j,'all');
% j_all_2(ii,jj) = sum(j2,'all');
% j_all_3(ii,jj) = sum(j3,'all');
% % j_all_1_1(ii) = sum(j1_1,'all');
% % j_all_2_1(ii) = sum(j2_1,'all');
% % j_all_3_1(ii) = sum(j3_1,'all');
% 
% I(ii,jj) =j_all_1(ii,jj)+j_all_2(ii,jj)+j_all_3(ii,jj);
% 
% 
% %I_1(ii) =j_all_1_1(ii)+j_all_2_1(ii)+j_all_3_1(ii);
% % 
%     end
% end
% 
% figure;
% surf(V_all,Vg_all,real(j_all_1)');
% figure;
% surf(V_all,Vg_all,real(j_all_2)');
% figure;
% surf(V_all,Vg_all,real(j_all_3)');
% 
% figure;
% plot(V_all,I(:,25))
% 
% % %%
% % %figure,plot(V_bias,I);
% % figure,plot(V_bias,movmean(I,15));
% % %figure,plot(V_bias,sgolayfilt(I,3,19));
% % 
% % 
% % 
% % %
% % % I_ft = (fft(I));
% % % zz = 1:1000;
% % % sigma = 50;
% % % gaussian = exp(-(zz/sigma).^2) + exp(-((1000-zz)/sigma).^2);
% % % figure,plot(V_bias,ifft(I_ft.*gaussian));
% % 
% % %
% % ax = gca;
% % ax.XAxisLocation = 'origin';
% % ax.YAxisLocation = 'origin';
% % 
% % [x0,y0] = ginput(1);
% % [minValue,closestIndex] = min(abs(V_bias-x0));
% % graphene_plot_band_structure(Vg,Vi(ind,closestIndex),deltaE_all(ind,closestIndex));
%% condutance


figure, imagesc(abs(diff(I_twist_new,1,1)'))
xlabel('Bias [V]')
ylabel('Gate voltage [V]')
set(gca,'YDir','normal')
figure,imagesc(abs(diff(I1(1:10:end,1:10:end),1,1)'))
xlabel('Bias [V]')
ylabel('Gate voltage [V]')
set(gca,'YDir','normal')


figure, plot(((I_twist_new(:,80))'))
%%

TER =0.5*( (Q_all2/(6*C0*er)).^2 -  (Q_all/(6*C0*er)).^2);
TER_as = repmat((-V_all*V_kp_0/12)',1,1003);
TER_as2 = repmat((-V_all*V_kp_0/18)',1,1003);

figure
surf(V_all,Vg_all,TER')
shading interp
hold on
%surf(V_all,Vg_all,TER_as')
surf(V_all,Vg_all,(TER_as2)')

shading interp
hold off
xlabel('V_b [V]')
ylabel('V_g [V]')
zlabel('TER')


figure,
plot(V_all,TER_as2(:,end))
hold on
plot(V_all,TER(:,end))
