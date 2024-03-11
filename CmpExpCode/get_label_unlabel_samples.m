function [X_label,gnd_label,X_unlabel,gnd_unlabel] = get_label_unlabel_samples(X,gnd,ind_label,ind_unlabel)
% return split of labeled and unlabeled samples according to the indices of
% labeled and unlabeled samples. The indices can be obtained by the
% functions 'get_label_unlabel_samples_indices.m' and 
% 'label_unlabel_split_index_producer.m'
n_view = length(X);

gnd_label = gnd(ind_label);
gnd_unlabel = gnd(ind_unlabel);
for v = 1:n_view
    X_label{v} = X{v}(ind_label,:);
    X_unlabel{v} = X{v}(ind_unlabel,:);
end
    
end
















