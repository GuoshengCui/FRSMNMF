function options = set_paras_dc_mvnmf_l1(dataset,options)
n_view = options.n_view;
X_label = options.X_label;
X_unlabel = options.X_unlabel;

options = choose_nn_dc_mvnmf_l1(dataset,options);
n_neighbors_lu = options.nn_lu;

for i = 1:n_view
    sigma_lu = 1;
    Si_lu = construct_W([X_label{i};X_unlabel{i}],n_neighbors_lu,sigma_lu);
    options.Slu{i} = Si_lu;
end

end

function [options] = choose_nn_dc_mvnmf_l1(dataset,options)

    if strcmp(dataset,'fei50')
        options.alpha = 1e5;
        options.beta = 1e4;
        options.nn_lu = 3;
    elseif strcmp(dataset,'orl40')
        options.alpha = 1e5;
        options.beta = 1e4;
        options.nn_lu = 2;
    else
        error('wrong dataset!!!');
    end
end


function W = construct_W(fea,num_knn,sigma)

      opts = [];
      opts.NeighborMode = 'KNN';
      opts.k = num_knn;
      opts.WeightMode = 'Binary';% HeatKernel Binary 
      opts.t = sigma;
      W = constructW(fea,opts);
      
end

