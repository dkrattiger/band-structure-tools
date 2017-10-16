function [V_damped,D_damped]=RPM(M,C,K,Vund,Dund,k1,k2)

%% Program that computes the damped eigenvalues and eigenvector from the
%% undamped ones via perturbation method

nmodes = size(Vund,2); %% Number of modes of interest
jmath=sqrt(-1);

%% Undamped System's Properties computed again this step can be avoided at
%% a later stage

%[Vund,Dund]=eig(K,M);
%% Sort the eigenvalues in ascending order

[D_sort,sort_indx]=sort(sqrt(diag(abs(Dund))));
V_sort=Vund(:,sort_indx);
D_und=diag(D_sort);


%% Scale the modes to unit modal mass

for ind=1:nmodes
    v=V_sort(:,ind);
    v=v/sqrt(v'*M*v);
    V_und(:,ind)=v;
end
%% Now apply the perturbation theory to predict Damped System's eigen
%% properties


%% Perturbation Method
omega =sqrt(diag(D_und));
Cp=V_und'*C*V_und; %% Damping matrix in modal co-ordinates


for n=1:nmodes
  D_damped(n,1)=D_und(n,n)+jmath*(Cp(n,n)/2);
%   D_damped(n)= ((-Cp(n,n) - sqrt(Cp(n,n)^2-4*D_und(n,n)))/2)^2;
    V_damped(:,n)=V_und(:,n);
    V_damped_imag(:,n)=zeros(size(V_damped(:,n)));
    for k=1:nmodes
        if(k~=n)
            alpha=jmath*Cp(k,n)*omega(n)/(omega(k)^2-omega(n)^2);
            V_damped_imag(:,n)= V_damped_imag(:,n)+alpha*V_und(:,k);
        end
    end
    V_damped(:,n)=V_und(:,n)+jmath*V_damped_imag(:,n);
    
end
