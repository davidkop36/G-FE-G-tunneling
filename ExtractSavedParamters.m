function [V_all,Vg_all,I1,Q_all,Q_bottom,Vi,qc,Vkp] =   ExtractSavedParamters(I_V)
V_all = I_V.V_all ;
Vg_all = I_V.Vg_all ;
I1 = I_V.I1;
Q_all = I_V.Q_all;
Q_bottom = I_V.Q_bottom;
Vi = I_V.Vi;
qc = I_V.qc;
Vkp = I_V.Vkp;
end