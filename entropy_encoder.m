function [ lossless_data ] = entropy_encoder( flat_data )
% ENTROPY_ENCODER funtion is based on an aritmethic encoder. It will take
% the fattened data and it will create an alphabeth, lossless_data, that
% will be send to the decoder. To know how many bits must be asign to a
% symbol, the probability of each siymbol must be calculated.
%
% INPUT VARIABLES:
%       'flat_data': Flattened residue after pre-processing
%
% OUTPUT VARIABLES:
%       'losslesss_data': Binari sequence encoded with n bits

n_frames = length(flat_data(1,:));
frame_length = length(flat_data(:,1));

% NOTE: We finally could not use the index value in the codification.
% Therefore, the probability template is fixed for symbols between 1 and 
% 5000.

%----------------------------- WITH INDEX (NOT USED) ----------------------
% %-------------------------------- MEAN ----------------------------------
% % For each frame of flattened prediction residue, the mean value of the 
% % frame, mu, is computed.
% mu = zeros(n_frames, 1);
% 
% for i = 1:n_frames
%     sum = 0;
%     for j = 1:frame_length
%         sum = sum + abs(flat_data(j,i));
%     end
%     mu(i) = sum / frame_length;
% end
% 
% %--------------------------- MEAN QUANTIZATION ----------------------------
% % The value mu is logarithmically quantized to an integer index.
% index = zeros(n_frames, 1);
% 
% for i = 1:n_frames
%     index(i) = floor(log2(mu(i)) + 0.5);
% end
% 
% % The first index is sent to the decoder in order to decode the first seq
% first_index = index(1);
% 
% %-------------------------- MEAN DE-QUANTIZATION --------------------------
% % The value index is locally dequantized as 2^index
% % 'deq_mu': dequanitzed mean of each frame
% deq_mu = zeros(n_frames, 1);
% 
% for i = 1:n_frames
%     deq_mu(i) = pow2(index(i));
% end
% 
% %--------------------------- ARITHMETIC ENCODER ---------------------------
% % 'lossless_data': binary arithmetic code
% % 'seq': Sequence of symbols array to be encoded
% 
% seq = zeros(frame_length+1, n_frames);
% lossless_data = cell(1,n_frames);
% 
% for i = 1:n_frames-1
%     % It is need to send the index value of future frame as the 
%     % first symbol of each sequence of symbols
%     seq(1,i) = index(i+1);
% 
%     % Flat_data would be the sequence of symbols to be encoded
%     for j = 1:frame_length
%         seq(j+1,i) = flat_data(j,i);
%     end
% 
%     s = min(min(seq(:,i))):1:max(max(seq(:,i)))+5;
%     
%     %------------------------------ SCALING -------------------------------
%     % The dequanitzed mean of the frame, which is used to scale a prob-
%     % ability template to generate a probability table for arithmetic coding
%     % 'scaled': Symbol scaled
%     
%     scaled = floor(s/deq_mu(i) + 0.5);
% 
%     % Avoid zero values by adding 1 due to arithenco only accept positive values
%     seq(:,i) = seq(:,i) + 1;
% 
%     %------------------------- PROBABILITY TEMPLATE -----------------------
%     % The probability template is a table containing a set of probability 
%     % density values, which are trained from a large amount of audio data.
%     prob_template = round(1e4*gaussmf(scaled,[0.3*max(max(scaled)) -0.1]));
%     
%     lossless_data{i} = arithenco(seq(:,i)', prob_template);
% end
% 
% % --------Last sequence-------- 
% % It is only sent the last sequence of symbols (without random index value)
% 
% % For example, the last index value equal to 1
% seq(1,n_frames) = 1;
% 
% % Flat_data would be the sequence of symbols to be encoded
% for j = 1:frame_length
%     seq(j+1,n_frames) = flat_data(j,n_frames);
% end
% 
% s = min(min(seq(:,n_frames))):1:max(max(seq(:,n_frames)))+5;
% 
% scaled = floor(s/deq_mu(n_frames) + 0.5);
% 
% % Avoid zero values by adding 1 due to arithenco only accept positive values
% seq(:,n_frames) = seq(:,n_frames) + 1;
% 
% prob_template = round(1e4*gaussmf(scaled,[0.3*max(max(scaled)) -0.1]));
%     
% lossless_data{n_frames} = arithenco(seq(:,n_frames)', prob_template);


%---------------------- WITHOUT INDEX (USED) ------------------------------
seq = zeros(frame_length, n_frames);
lossless_data = cell(1,n_frames);

%------------------------- PROBABILITY TEMPLATE -----------------------
% The probability template is a table containing a set of probability 
% density values, which are trained from a large amount of audio data.
s = 1:1:5000;
prob_template = round(1e4*gaussmf(s,[0.3*max(max(s)) -0.1]));

for i = 1:n_frames
    % Flat_data would be the sequence of symbols to be encoded
    for j = 1:frame_length
        seq(j,i) = flat_data(j,i);
    end
    
    % Avoid zero values by adding 1 due to arithenco only accept positive values
    seq(:,i) = seq(:,i) + 1;
    
    lossless_data{i} = arithenco(seq(:,i)', prob_template);
end

end

