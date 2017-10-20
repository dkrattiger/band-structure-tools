function [kappa_sorted,i_sort] = curve_sort_k_w(kappa)

[n_curves,n_om] = size(kappa);

kappa_sorted = zeros(size(kappa));
i_sort = zeros(size(kappa));
% kappa_slice = kappa(:,1);


kappa_sorted(:,1) = kappa(:,1);
i_sort(:,1) = 1:n_curves;

% figures
for i = 2:n_om
    
    % get previous kappa_slice
    kappa_slice_prev = kappa_sorted(:,i-1);
%     [~,i_sort_mag_last] = sort(kappa_slice_prev);
    
    kappa_slice = kappa(:,i);
    i_sort_slice = zeros(n_curves,1);
    
    i_left = 1:length(kappa_slice);
    i_left_prev = 1:length(kappa_slice_prev);
    
    [kappa_slice_prev_grid,kappa_slice_grid] = meshgrid(kappa_slice_prev,kappa_slice);
    diff_mat = kappa_slice_prev_grid-kappa_slice_grid;
    for j = 1:n_curves
%         diff_vec = abs(kappa_slice_prev(i_sort_mag_last(j))-kappa_slice);        
        [diff_mat_mins,i_mins] = min(diff_mat);
        [~,j_min] = min(diff_mat_mins);
        i_sort_slice(i_left_prev(j_min)) = i_left(i_mins(j_min));
        
        diff_mat(i_mins(j_min),:) = [];
        diff_mat(:,j_min) = [];
        
        i_left(i_mins(j_min)) = [];
        i_left_prev(j_min) = [];
    end
    
    kappa_sorted(:,i) = kappa(i_sort_slice,i);
    i_sort(:,i) = i_sort_slice;
end

% for i = 2:n_om
%     
%     % get previous kappa_slice
%     kappa_slice_prev = kappa_sorted(:,i-1);
%     [~,i_sort_mag_last] = sort(kappa_slice_prev);
%     
%     kappa_slice = kappa(:,i);
%     i_sort_slice = zeros(n_curves,1);
%     i_left = 1:length(kappa_slice);
%     for j = 1:n_curves
%         diff_vec = abs(kappa_slice_prev(i_sort_mag_last(j))-kappa_slice);        
%         [~,i_min] = min(diff_vec);
%         i_sort_slice(i_sort_mag_last(j)) = i_left(i_min);
%         kappa_slice(i_min) = [];
%         i_left(i_min) = [];
%     end
%     
%     kappa_sorted(:,i) = kappa(i_sort_slice,i);
%     i_sort(:,i) = i_sort_slice;
%     
% %     plot3(real(kappa_sorted(:,1:i)),imag(kappa_sorted(:,1:i)),1:i,'.-')
% %     pause
% end