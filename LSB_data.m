function [ LSB_data ] = LSB_data(q_parq, pred_residues)
%LSB_DATS Summary of this function goes here
%   Detailed explanation goes here

n_residues = length(q_parq(1,:));
PARCOR_order=length(q_parq(:,1));
%residue_length=length(pred_residues(:,1));

n_shifts=zeros(PARCOR_order,n_residues);
n_shift=zeros(PARCOR_order,n_residues);
%flat_data=zeros(residue_length,n_residues);
%LSB_data=zeros(PARCOR_order,n_residues);

for k=1:n_residues
  sum=0;  
    for i=0:PARCOR_order-1
    
           for j=i+1:PARCOR_order
        
                
                n_shift(j,k) = log2(1/(1-q_parq(j,k)^2))+0.5+sum;
                n_shifts=floor(n_shift);
                sum=n_shift(j,k);
            end
    end
          
end

LSB_data=mod(pred_residues,2.^n_shifts);

end

