clear all

%% Input parameters of the model
er = 3; %FE material dielectric constant
er_gate = 3.8; %dielectric constant of the gate material
d = 5e-9; % Junction thickness. Default to 5nm . Try 1nm, per Moshe's recommendation 
d_gate = 90e-9; %gate thickness, 15e-9 standard 
theta  =0; % 2 *pi/180; %twist angle in between two graphene sheets
V_kp_0 = +0.1; % Bare ferroelctricity voltage
qc = (12e-9)^(-1); % Scattering potential shape paramaterer (value according to Brittnel)

% Bias_gate range of interest
V_b_max = 2;
V_b_min = -2;
V_g_max = 30;
V_g_min = -30;

course_grid = [101, 103]; %used for electrostatics calculation
fine_grid = [1001, 1003]; %used for I-V current calculation

% plot customization
plot_type_option = 1; %set to 1 to plot I_{+} on the color plot or to 2 to plot I_{+} + I_{-}

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

% bottom sheet charges

Q_bottom= Q1_all-Q_all;
Q_bottom_2 = Q1_all2-Q_all2;
Q_bottom_3 = Q1_all3-Q_all3;

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



%% Asymptotic expression vs the numeric results
indg =floor(0.6*size(V_all,2));
Q0 = (V_all+V_kp_0)*er*C0;
Q0_prim = Vg_all(indg)*C1;

f_fit2 =  beta*sqrt(er*C0./V_all);

figure;
f_thing_zero = f_thing(:,indg);
plot(V_all,f_thing(:,indg),'LineWidth',2);
xlabel('V_b [V]','FontSize',20)
ylabel(' f','FontSize',20)
title(['For V_g = ',num2str(Vg_all(indg),'%.2f'),' V'  ],'FontSize',20)
hold on
plot(V_all((V_all)>V_kp_0*4),f_fit2((V_all)>V_kp_0*4),'r--','LineWidth',2)

f_fit_peak = 2./(sqrt(1+4*abs(V_all+V_kp_0)./(C0*er*beta^2)) +sqrt(1+4*abs(V_all-V_kp_0)./(C0*er*beta^2))  );


%% Find Fermi energies relative to Dirac point
% Polarization up
Ef_top = -sqrt(4*pi*abs(Q_all)/e_el).*sign(Q_all)*hbar*Vf/(2*e_el);
Ef_bott = -sqrt(4*pi*abs(Q_bottom)/e_el).*sign(Q_bottom)*hbar*Vf/(2*e_el);
V_all_rep = repmat(V_all',1,length(Vg_all));
Vi = V_all_rep + Ef_top - Ef_bott;


gamma = e_el*(Ef_bott-Ef_top-V_all_rep)/(hbar*Vf*qc);
kt1 = Ef_top*e_el/(hbar*Vf*qc);
kt2 = (Ef_top+V_all_rep)*e_el/(hbar*Vf*qc);
k_vi = e_el*Vi/(hbar*Vf*qc);

% Polarization down 
Ef_top_2 =- sqrt(4*pi*abs(Q_all2)/e_el).*sign(Q_all2)*hbar*Vf/(2*e_el);
Ef_bott_2 = -sqrt(4*pi*abs(Q_bottom_2)/e_el).*sign(Q_bottom_2)*hbar*Vf/(2*e_el);
Ef_top_3 =- sqrt(4*pi*abs(Q_all3)/e_el).*sign(Q_all3)*hbar*Vf/(2*e_el);
Ef_bott_3 =- sqrt(4*pi*abs(Q_bottom_3)/e_el).*sign(Q_bottom_3)*hbar*Vf/(2*e_el);
Vi_2 = V_all_rep + Ef_top_2 - Ef_bott_2;
Vi_3 = V_all_rep + Ef_top_3 - Ef_bott_3;
gamma_2 = e_el*(Ef_bott_2-Ef_top_2-V_all_rep)/(hbar*Vf*qc);
kt1_2 = Ef_top_2*e_el/(hbar*Vf*qc);
kt2_2 = (Ef_top_2+V_all_rep)*e_el/(hbar*Vf*qc);

k_vi_2 = e_el*Vi_2/(hbar*Vf*qc);

reduced_twist_vector = (8*pi/(3*1.42e-10))*sin(theta/2)/qc;
%% current calculation
I1 = Calculate_current_graphene(kt1,kt2,k_vi,reduced_twist_vector);
I2 = Calculate_current_graphene(kt1_2,kt2_2,k_vi_2,reduced_twist_vector);
%%
if theta == 0
gde = floor(0.75*size(V_all,2));
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
end

%% Plot the two I-V curves for a given gate on top of each other
gde1 =floor(0.55*size(V_all,2));

figure;
plot(V_all,I1(:,gde1),'LineWidth',2)
hold on
plot(V_all,I2(:,gde1),'LineWidth',2)
xlabel('Bias [V]','FontSize',20)
ylabel('Current [a.u.]','FontSize',20)
title(['Current vs. bias for Vg =  ',num2str(Vg_all(gde1),'%.2f'),'V'],'FontSize',20 )
legend('V_{KP}^{(0)} > 0 ','V_{KP}^{(0)} < 0 ' ,'Location','bestoutside','FontSize',15)
ax = gca;
ax.FontSize = 15;

% exatract lines on 2D graphs corresponding to physically signficant cases
Vi_zero_locus = ExtractMinimumInDirection(abs(Vi-reduced_twist_vector*qc*hbar*Vf/e_el),1);
Vi_zero_locus_other = ExtractMinimumInDirection(abs(Vi+reduced_twist_vector*qc*hbar*Vf/e_el),1);
Qtop_zero_locus = ExtractMinimumInDirection(abs(Q_all),1);
Qbottom_zero_locus = ExtractMinimumInDirection(abs(Q_bottom),1);
Vi_zero_locus_2 = ExtractMinimumInDirection(abs(Vi_2-reduced_twist_vector*qc*hbar*Vf/e_el),1);
Vi_zero_locus_other_2= ExtractMinimumInDirection(abs(Vi_2+reduced_twist_vector*qc*hbar*Vf/e_el),1);
Qtop_zero_locus_2 = ExtractMinimumInDirection(abs(Q_all2),1);
Qbottom_zero_locus_2 = ExtractMinimumInDirection(abs(Q_bottom_2),1);
Vi_zero_locus_0 = ExtractMinimumInDirection(abs(Vi_3-reduced_twist_vector*qc*hbar*Vf/e_el),1);
Qtop_zero_locus_0 = ExtractMinimumInDirection(abs(Q_all3),1);
Qbottom_zero_locus_0 = ExtractMinimumInDirection(abs(Q_bottom_3),1);
Qprim_zero_locus_0 = ExtractMinimumInDirection(abs(Q1_all),1);


% Qcompens_zero_locus_0 = ExtractMinimumInDirection(abs(Q_all + Q_all2),1);

% For validation purposes one can compare the theoretical and numerical
% curves
[V_bias_zero_Q_bottom,V_gate_zero_Q_top,V_gate_zero_Vi]= CalculatePlotTheoreticalCurves(V_all,Vg_all,V_kp_0,C0,C1,beta,er,Qbottom_zero_locus,Qtop_zero_locus,Vi_zero_locus,theta);
[V_bias_zero_Q_bottom2,V_gate_zero_Q_top2,V_gate_zero_Vi2]= CalculatePlotTheoreticalCurves(V_all,Vg_all,-V_kp_0,C0,C1,beta,er,Qbottom_zero_locus_2,Qtop_zero_locus_2,Vi_zero_locus_2,theta);
%[V_bias_zero_Q_bottom0,V_gate_zero_Q_top0,V_gate_zero_Vi0]= CalculatePlotTheoreticalCurves(V_all,Vg_all,0,C0,C1,beta,er,Qbottom_zero_locus_0,Qtop_zero_locus_0,Vi_zero_locus_0);


%%



figure;
if plot_type_option==2
    %plot the total equidomenial current in some range
    imagesc(V_all,Vg_all,(I2')+(I1'));

else
    % plot just the positive orientation current onto 2D plot
    imagesc(V_all,Vg_all,(I1'));
end
set(gca,'YDir','normal')

hold on
lw_def = 2;

%theoretical curves plot
plot(V_all,V_gate_zero_Q_top,'g--','LineWidth',lw_def)
plot(V_bias_zero_Q_bottom,Vg_all,'y--','LineWidth',lw_def)
if theta ==0
    %plot(V_all,Vi_zero_locus,'r--','LineWidth',lw_def)
      plot(V_all(Vi_zero_locus),Vg_all,'r--','LineWidth',lw_def)
else
    plot(V_all(Vi_zero_locus),Vg_all,'r--','LineWidth',lw_def)
    plot(V_all(Vi_zero_locus_other),Vg_all,'r--','LineWidth',lw_def)
end    

plot(V_all,V_gate_zero_Q_top2,'g.','LineWidth',lw_def)
plot(V_bias_zero_Q_bottom2,Vg_all,'y.','LineWidth',lw_def)
if theta==0
    %plot(V_all,Vi_zero_locus_2,'r--','LineWidth',lw_def)
    plot(V_all(Vi_zero_locus_2),Vg_all,'r.','LineWidth',lw_def)
else
    plot(V_all(Vi_zero_locus_2),Vg_all,'r.','LineWidth',lw_def)
    plot(V_all(Vi_zero_locus_other_2),Vg_all,'r.','LineWidth',lw_def)

end 

xlabel('Bias [V]','FontSize', 20)
ylabel('Gate [V]','FontSize', 20)
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
% Pick points where to a make a 1D cut
closestIndex2 = floor(0.834*size(Vg_all,2)); %Fix at 20V gate at 1003 points
midindex = floor(0.66*size(Vg_all,2)); %Fix at 9.52V gate at 1003 points
yline(Vg_all(closestIndex2));
yline(Vg_all(midindex));

if theta == 0
%select points to plot the bandstructure at, default 6
[x0,~] = ginput(6); 
scatter(x0,repmat(Vg_all(closestIndex2),1,6))
scatter(x0,repmat(Vg_all(midindex),1,6))
else
%select points to plot the bandstructure at, default 4
[x0,~] = ginput(4);
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
% 2nd line is the midindex plot
else
graphene_plot_band_structure_qc_twisted_debug(V_all(closestIndex),Vg_all(closestIndex2),Vi(closestIndex,closestIndex2),Q_all(closestIndex,closestIndex2),Q_bottom(closestIndex,closestIndex2),theta);
end
end

%% Plot the I-V curves at those gates
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


%% Distance between the peaks
if theta ==0
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

end
