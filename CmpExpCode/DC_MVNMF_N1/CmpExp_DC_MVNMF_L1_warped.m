
function [ac_m,ac_s,nmi_m,nmi_s,ar_m,ar_s] = ...
    CmpExp_DC_MVNMF_L1_warped(options)

n_view = options.n_view;
num_nmf = options.num_nmf;% num conduct nmf
num_km = options.num_km;% num repeat of kmeans
nClas = options.nClas;
gnd_perm = options.gnd_perm;
options.gnd_known = options.gnd_label;

for v = 1:n_view
    fea{v} = [options.X_label{v};options.X_unlabel{v}];
end

AC = [];NMI = [];AR = [];
Fscore = [];Precision = [];Recall = [];
rand('seed',2);
for inmf = 1:num_nmf

    tic
    [W,H,Hc,dc_mvnmf_L1_obj] = DC_MVNMF_L1_Norm(fea,options);
    toc
    %% 1/3*(H{1}+H{2}+H{3}) used!!!
    Ht = 0;
    for i = 1:options.n_view
        Ht = Ht + H{i};
    end
    Ht = 1/n_view*Ht;

    label = kmeans(Ht' ,nClas, 'Replicates',num_km,'Distance','correlation');
    [ac,nmi,~] = result(label,gnd_perm);
    [ar,~,~,~] = RandIndex(gnd_perm,label);
    [f,p,r] = compute_f(gnd_perm,label);
    AC = [AC,ac];
    NMI = [NMI,nmi];
    AR = [AR,ar];
    Fscore = [Fscore, f];
    Precision = [Precision, p];
    Recall = [Recall, r];
    disp(strcat('(',num2str(roundn(ac,-3)),',',num2str(roundn(nmi,-3)),...
        ',',num2str(roundn(ar,-3)),',',num2str(roundn(f,-3)),...
        ',',num2str(roundn(p,-3)),',',num2str(roundn(r,-3)),')'));
end
ac_m = mean(AC);
ac_s = std2(AC);
nmi_m = mean(NMI);
nmi_s = std2(NMI);
ar_m = mean(AR);
ar_s = std2(AR);
fprintf('dc_mvnmf_l1 the accuracy is : %4.2f+%2.1f\n',mean(AC)*100,std2(AC)*100);
fprintf('dc_mvnmf_l1 the NMI is : %4.2f+%2.1f\n',mean(NMI)*100,std2(NMI)*100);
fprintf('dc_mvnmf_l1 the AR is : %4.2f+%2.1f\n',mean(AR)*100,std2(AR)*100);
% fprintf('dc_mvnmf_l1 the Fscore is : %4.2f+%2.1f\n',mean(Fscore)*100,std2(Fscore)*100);
% fprintf('dc_mvnmf_l1 the Precision is : %4.2f+%2.1f\n',mean(Precision)*100,std2(Precision)*100);
% fprintf('dc_mvnmf_l1 the Recall is : %4.2f+%2.1f\n',mean(Recall)*100,std2(Recall)*100);

