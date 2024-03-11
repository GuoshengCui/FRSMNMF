function [AC,MI,label] = result(label,gnd)

label = bestMap(gnd,label);
AC = length(find(gnd == label))/length(gnd);
MI = MutualInfo_cai(gnd,label);
% MI = MutualInfo_yang(gnd,label);
% MI = NMI(gnd,label);
% [AC,MI] = AcMI(label,gnd);

