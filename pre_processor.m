function [ flat_data, LSB_data] = pre_processor( q_parq, residues )
% Pre_processor function takes the prediction residues and it apply them 
% a down-shift by a number of BITS, so their amplitude decrease and the 
% envelop of the prediction resiues is flattened.
%
% INPUT VARIABLES:
%       'q_parq': Quantized PARCOR coefficients
%                 Dimension Matrix: (parcor_order,n_frames)
%                                   (i.e. parcor_order = 2^7bits = 128 )
%       'residues': Prediction residue
%
% OUTPUT VARIABLES:
%       'flat_dat': Flattened residue after pre-processing
%       'LSB_data': Bits that are removed when the down-shift operation
%                   is done

%--------------------------------------------------------------------------

%   'n_frames': Number of frames processed
%   'residue_length': Residue length in one frame

n_frames = length(residues(1,:));
parcor_order = length(q_parq(:,1));
frame_length = length(residues(:,1));

%---------------------- NUMBER OF SHIFTS-----------------------------------
% 'shift': Number of down-shifts applied to each residue.
%           Shifts number are calculated from the quantized PARCOR
%           coefficients (Bits that we have to decrease from the
%           prediction residues from each frame).

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

%---------------------------- DOWN-SHIFT OPERATION ------------------------

          %--------------- LSB DATA -----------------
% LSB_data contains information about the sign of the signal    

LSB_data = zeros(L,n_frames);

for i = 1:n_frames
    for j = 1:L
        LSB_data(j,i) = mod(abs(residues(j,i)), pow2(shift(j,i)));
    end
end
          %------------ FLATTENED RESIDUES ----------
% Flattened_residues contains information about the sign of the signal 
% although the standard does not. That is because some floor operations in 
% LSB_data computing looses the sign due to the value is 0.

flat_data = zeros(frame_length,n_frames);

for i = 1:n_frames
    for j = 1:L
        dec = abs(residues(j,i));
        
        % bitshift(A,k) returns A shifted to the left by k bits, equivalent 
        % to multiplying by 2k. Negative values of k correspond to shifting 
        % bits right or dividing by 2|k| and rounding to the nearest integer 
        % towards negative infinity.
        
        flat_data(j,i) = bitshift(dec, -shift(j,i));
    end

    for j = L+1:frame_length
        flat_data(j,i) = abs(residues(j,i));
    end
end

end

