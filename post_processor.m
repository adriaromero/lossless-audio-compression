function [ post_residues ] = post_processor( signs, dflat_data, LSB_data, q_parq )
% With the post_processor function we reconstruct the prediction residues
% using the flattened resitues that have been encoded and decoded using
% entropy coding and the LSB_data that have been transmitted directly from
% the pre-processing. Therefore an Up-shift operation must be done and to
% compute the number of shifts we use the PARCOR quantized coefficients. 
%
% INPUT VARIABLES:
%
%       'dflat_data': Decoded flattened residue 
%       'LSB_data': Bits that are add when the up-shift operation
%                   is done.
%       'q_parq': Quantized PARCOR coefficients
%
% OUTPUT VARIABLES:
%
%       'post_residues': Post processing prediction residues

%--------------------------------------------------------------------------
%   'n_frames': Number of frames processed

n_frames = length(dflat_data(1,:));
parcor_order = length(q_parq(:,1));
frame_length = length(dflat_data(:,1));

%---------------------- NUMBER OF SHIFTS ----------------------------------
% 'shift': Number of down-shifts applied to each residue.
%          Shifts number are calculated from the quantized PARCOR
%          coefficients (Bits that we have to decrease from the
%          prediction residues from each frame).

L = parcor_order;
shift = zeros(L,n_frames);
sum1 = 0;
sum2 = 0;
RA_shift12 = dlmread('RA_shift12.txt');   % RA_shift12 table
RA_shift = dlmread('RA_shift.txt');   % RA_shift table

for i = 1:n_frames
    for j = 1:2
        for k = 1:j
            sum1 = sum1 + RA_shift12(RA_shift12(:,1) == q_parq(k,i),2);
        end
        shift(j,i) = floor((2.^12 + sum1)/2.^13);
        sum1 = 0;
    end

    for j = 3:L
        for k = 1:2
            sum1 = sum1 + RA_shift12(RA_shift12(:,1) == q_parq(k,i),2);
        end
        
        for k = 3:j
            sum2 = sum2 + RA_shift(RA_shift(:,1) == abs(q_parq(k,i)),2);
        end
        
        shift(j,i) = floor((2.^12 + sum1 + sum2)/2.^13);
        sum1 = 0;
        sum2 = 0;
    end
end

%---------------------- POST PREDICTION RESIDUE ---------------------------
post_residues = zeros(frame_length, n_frames);

for i = 1:n_frames
    for j = 1:L
        post_residues(j,i) = dflat_data(j,i) * pow2(shift(j,i)) + LSB_data(j,i);
        
        if (signs(j,i) == -1)   % Correction of sign value
            post_residues(j,i) = - post_residues(j,i);
        end
    end
       
    for j = L+1:frame_length
        post_residues(j,i) = dflat_data(j,i);
        
        %Using entropy coders
        if (signs(j,i) == -1)   % Correction of sign value
            post_residues(j,i) = - post_residues(j,i);
        end
    end 
end

end

