function [Yeti, Zeti] = electrode_skin_interface(f)

%Data extracted from S. Kaufmann et al., "Measurements of Electrode Skin
%Impedance using Carbon Rubber Electrodes - First Results", Journal of
%Physics: Conference Series, 434: 1-4, 2013. 

%%Dry-skin, pre-gelled electrodes. Fig 2-e,f
Z = [76
70
64
58
52
48
44
41
39
37
35
33
32
30
29
27
26
26];  %Ohms

Z_ang = [-55
-54
-51
-48
-44
-42
-39
-36
-34
-32
-30
-28
-25
-24
-21
-18
-16
-14
];

freq = [12
14
16
20
25
29
35
42
49
59
69
81
98
115
148
196
244
296]*1e3; %Hz
R_measured = 20e-3;
Zeti = interp1(freq,abs(Z).*exp(1i*deg2rad(Z_ang)),f,'linear','extrap');
Yeti = 1./Zeti;

Yeti = Yeti/(pi*R_measured^2); %normalized per area