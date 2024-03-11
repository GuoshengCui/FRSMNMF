function options = set_paras_common(options,path,dataset,tail_fix,idata)


    

load(strcat(path,'\',dataset,'_label_unlabel_ind_split_',tail_fix,'_',...
    num2str(idata),'.mat'));
[X_label,gnd_label,X_unlabel,gnd_unlabel] = ...
    get_label_unlabel_samples(options.fea,options.gnd,ind_label,ind_unlabel);
options.X_label = X_label;
options.X_unlabel = X_unlabel;
options.gnd_label = gnd_label;
options.gnd_known = gnd_label;
options.gnd_unlabel = gnd_unlabel;
options.gnd_perm = [gnd_label;gnd_unlabel];
options.nClas = length(unique(options.gnd));
options.n_view = length(options.fea);

