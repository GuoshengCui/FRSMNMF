
close all
clear
clc
% dataset = {'orl40','fei50'};
dataset = {'fei50'};
ratio = 0.3;
nSubSpace = 1;
options.ratio = ratio;
%     options.num_km = 10;
%     options.num_nmf = 10;
options.max_iter = 300;
options.num_km = 1;
options.num_nmf = 1;
%     options.max_iter = 1;
options.nSubSpace = 1;

for idataname = 1:length(dataset)
    [fea,gnd] = ld_image_mv_dataset(dataset{idataname});
    options.gnd = gnd;
    options.fea = fea;
    path = strcat(pwd,'\Loadata\semi-supervised-data-split\',dataset{idataname},'\');
    tail_fix = choose_tail_fix(ratio);
    
    for idata = 1:5
        % set common paras
        options = set_paras_common(options,path,dataset{idataname},tail_fix,idata);
        % set paras of DC_MVNMF_L2.
        options = set_paras_dc_mvnmf_l2(dataset{idataname},options); % 'lu' denotes 'label' and 'unlabel'.
        
        % DC_MVNMF_L2
        [ac_m,ac_s,nmi_m,nmi_s,ar_m,ar_s] = CmpExp_DC_MVNMF_L2_warped(options);
        AC(idata) = ac_m;
        NMI(idata) = nmi_m;
        AR(idata) = ar_m;
    end
    rest_DC_MVNMF_L2{1} = [mean(AC)*100,std2(AC)*100];
    rest_DC_MVNMF_L2{2} = [mean(NMI)*100,std2(NMI)*100];
    rest_DC_MVNMF_L2{3} = [mean(AR)*100,std2(AR)*100];
    fprintf('DC_MVNMF_L2 the accuracy is : %4.2f+%2.1f\n',mean(AC)*100,std2(AC)*100);
    fprintf('DC_MVNMF_L2 the NMI is : %4.2f+%2.1f\n',mean(NMI)*100,std2(NMI)*100);
    fprintf('DC_MVNMF_L2 the AR is : %4.2f+%2.1f\n',mean(AR)*100,std2(AR)*100);
end

function tail_fix = choose_tail_fix(ratio)

if ratio==0.1
    tail_fix = 'r1_10';
elseif ratio==0.2
    tail_fix = 'r2_10';
elseif ratio == 0.3
    tail_fix = 'r3_10';
end
end
