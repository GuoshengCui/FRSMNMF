function [fea,gnd] = ld_image_mv_dataset(dataset)
% output fea: cell type {n*m1,n*m2,...,n*mk}
% output gnd: n*1


%% load dataset 
switch lower(dataset)
    case 'orl40'
        load(strcat(pwd,'\mv_orl_40.mat'));
        X{1} = pixel;X{2} = gabr;X{3} = lbp;
        n_view = length(X);
        nClas = length(unique(gnd));
        fea = X;
    case 'fei50' 
        load(strcat(pwd,'\mv_fei_50.mat'));
        X{1} = pixel;X{2} = gabr;X{3} = lbp;
        n_view = length(X);
        nClas = length(unique(gnd));
        fea = X;
    otherwise
        error('wrong dataset!');
end
%% normalize dataset 
n_view = length(fea);
fea = normalize_fea(fea,n_view,0,1);% maxmin:1;else:norm sample to unit vec
%% disp info of dataset 
disp(' ');
disp(strcat('''',dataset,'''',' dataset info: '));
disp('//------------------------------------------------//');
disp('nClas: ');
nClas = length(unique(gnd));
disp(strcat('-----',num2str(nClas)));
for i = 1:n_view
    disp(strcat('nSmp x mFea: '));
    disp(strcat('-----(',num2str(size(fea{i},1)),', ',...
        num2str(size(fea{i},2)),')'));
end
str = [];
for i = 1:length(unique(gnd))
    if i == 1
        str = strcat('-----(',num2str(i),':',num2str(sum(gnd==i)),')');
    else
        str = strcat(str,'---','(',num2str(i),':',num2str(sum(gnd==i)),')');
    end
end
disp('sample of each class: ');
disp(str);
disp('//------------------------------------------------//');
end
%% func to normalize dataset
function [X] = normalize_fea(X,n_view,maxmin,unit_vector)
    % if maxmin=1, unit_vector has no influence.
    if maxmin
        for i = 1:n_view
            X_temp = full(X{i});
            max_col = max(X_temp,[],1);
            min_col = min(X_temp,[],1);
            range_col = max_col - min_col;
            % eliminate dummy features
            ind_dummy_f = find(range_col==0);
            X_temp(:,ind_dummy_f) = [];
            min_col(ind_dummy_f) = [];
            range_col(ind_dummy_f) = [];
            % 
            diff_mat = X_temp - min_col;
            X{i} = diff_mat./range_col;
        end
    else
        for i = 1:n_view
            % unit vector
            X_temp_2 = X{i};
            if unit_vector
                norm_col = sqrt(sum(X_temp_2.^2,2));
            else
                norm_col = sum(X_temp_2,2);
            end
            % 
            X{i} = X_temp_2./norm_col;
        end
    end
    
end
