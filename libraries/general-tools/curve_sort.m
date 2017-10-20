function [w_sort,PHI_sort,i_sort] = curve_sort(w,PHI);


% array dimensions
[n_dof,n_curves,n_kap] = size(PHI);

% % normalize mode shapes to unit length
% vec_lengths = ones(n_dof,1)*sqrt(sum(PHI.*conj(PHI),1));
% PHI = PHI./vec_lengths;

% preallocate sorted arrays
i_sort = zeros(n_curves,n_kap);
PHI_sort = zeros(n_dof,n_curves,n_kap);
w_sort = zeros(n_curves,n_kap);

% don't change sort at first k-point
i_sort(:,1) = 1:n_curves;
PHI_sort(:,:,1) = PHI(:,:,1);
w_sort(:,1) = w(:,1);

for i = 2:n_kap
    
    % normalize mode vectors to unit magnitude (this step is inefficient)
    PHI1 = PHI_sort(:,:,i-1)*diag(diag(PHI_sort(:,:,i-1)'*PHI_sort(:,:,i-1)).^(-0.5));
    PHI2 = PHI(:,:,i)*diag(diag(PHI(:,:,i)'*PHI(:,:,i)).^(-0.5));
    
    % orthonormalize eigenvector sets
    [PHI1,~] = qr(PHI1,0);
    [PHI2,~] = qr(PHI2,0);
    
    % compute modal assurance criterion
    MAC = abs(PHI1'*PHI2);
    
    for j = 1:n_curves
        
        % find maximum entry in every column
        [maxes,i_max_rows] = max(MAC);
        
        % find column with maximum entry
        [~,i_max_col] = max(maxes);
        
        % row with maximum entry
        i_max_row = i_max_rows(i_max_col);
        
        % place maximum column index into maximum row index entry of i_sort
        i_sort(i_max_row,i) = i_max_col;
        
       
        % remove current row and column from consideration
        MAC(i_max_row,:) = nan;
        MAC(:,i_max_col) = nan;
        
    end

    % re-sort PHI and w to match previous set
    PHI_sort(:,:,i) = PHI2(:,i_sort(:,i));
    w_sort(:,i) = w(i_sort(:,i),i);
    
    % align phase of mode shapes
    for j = 1:n_curves
        PHI_sort(:,j,i) = PHI_sort(:,j,i)/(PHI_sort(:,j,i-1)'*PHI_sort(:,j,i));
    end
end


    