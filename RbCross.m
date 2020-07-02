function sig_rb = RbCross(v_ion)
% v_ion in cm/s
% sig_rb in cm^2

a = 42E-8;
b = 1.85E-8;

sig_rb = (a-b*log(v_ion)).^2;

%sig_rb = 9.0441e-13;