function [ audio_output ] = reconstructor( q_parq, post_residues)
% RECONSTRUCTOR 
% In audio_ouput, the quantized PARCOR coefficients are extracted from 
% the bitstream, dequantized, and converted to linear predictor coefficients, 
% which are identical to those used in the encoder. The linear predictor 
% generates a prediction, which is added to the decoded prediction residue 
% to reconstruct the original input sample.
%
% INPUT VARIABLES:
%       'post_residues': Prediction residue
%       'q_parq':   Quantized PARCOR coefficients. 
%                   Matrix dimensions: (parcor_order,n_frames)
%
% OUTPUT VARIABLES:
%       'reconstruction': Audio output signal with format .wav, .ogg, .flac,
%                       .mp3 or .mp4
%

%--------------------------------------------------------------------------
%   'n_frames': Number of frames processed

n_frames = length(q_parq(1,:));
parcor_order = length(q_parq(:,1));
frame_length = length(post_residues(:,1));
n_samples = n_frames * frame_length;

%----------------------------PARCOR DE-QUANTIZATION------------------------
% 'deq_parq': De-quantizated PARCOR coefficients

deq_parq = zeros(parcor_order, n_frames);

for i = 1:n_frames
    deq_parq(1,i) = (2 * ((exp(q_parq(1,i)/64*log(3/2))-(2/3)) * 6/5).^2)-1; 
    deq_parq(2,i) = -(2 * ((exp(q_parq(2,i)/64*log(3/2))-(2/3)) * 6/5).^2)+1;
    deq_parq(3:parcor_order,i) = q_parq(3:parcor_order,i)/64; 
end

%----------------------------PARCOR TO LPC---------------------------------
% The reconstructed PARCOR coefficients are converted to k-order 
% (1 < j < lpc_order) LPC coefficients lpc(j,1..j) 
% 'lpc' : Linear Prediction Coefficients
lpc = zeros(parcor_order, parcor_order, n_frames);
 
for j = 1:n_frames
       for m = 1:parcor_order
            lpc(m,m,j) = deq_parq(m,j);
   
            for i = 1:m-1
                lpc(i,m,j) = lpc(i,m-1,j) + deq_parq(m,j)*lpc(m-i,m-1,j);
            end
       end
end

%-------------------------LINEAR PREDICTOR DECODER-------------------------
% The Linear Predictor Decoder reconstructs the audio input signal from the
% prediction post_residues of each sample in the frame.
% 'audio_ouput': audio_ouput of prediction post_residues

reconstruction = zeros(frame_length, n_frames);
% 
% for i = 1:n_frames
%     reconstruction(1,i) = post_residues(1,i);
%     
%     for j = 2:parcor_order
%         sum = 0;
%         for k = 1:j-1
%             sum = sum + lpc(k,j,i)*reconstruction(j-k,i);
%         end 
%         reconstruction(j,i) = post_residues(j,i)-round(sum);
%     end
%    
%     for j = parcor_order+1:frame_length
%         sum = 0;
%         for k = 1:parcor_order
%             sum = sum + lpc(k,parcor_order,i)*reconstruction(j-k,i);
%         end 
%         reconstruction(j,i) = post_residues(j,i)-round(sum);
%     end
% end

for i = 1:n_frames
    reconstruction(1,i) = post_residues(1,i);
    
    for j = 2:parcor_order
        sum = 0;
        for k = 1:j-1
            sum = sum + lpc(k,j,i)*reconstruction((j-k),i);
        end 
        reconstruction(j,i) = post_residues(j,i)-round(sum);
    end
   
    for j = parcor_order+1:frame_length
        sum = 0;
        for k = 1:parcor_order
            sum = sum + lpc(k,parcor_order,i)*reconstruction(j-k,i);
        end 
        reconstruction(j,i) = post_residues(j,i)-round(sum);
    end
end

audio_output = reshape(reconstruction/2.^15, 1, 400*176);

end

