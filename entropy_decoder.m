function [ dflat_data ] = entropy_decoder( lossless_data )
% ENTROPY_DECODER funtion is based on an aritmethic decoder. It will take
% the lossless data and it will decode the binary arithmetic in the 
% vector code to recover the corresponding sequence of flat_data.
%
% INPUT VARIABLES:
%       'first_index': Index of the first sequence to be decoded
%       'losslesss_data': Binari sequence encoded with n bits
%
% OUTPUT VARIABLES:
%       'dflat_data': Decodified Flattened residue after pre-processing

n_frames = length(lossless_data(1,:));
frame_length = 400;

% NOTE: We finally could not use the index value in the decodification.
% Therefore, the probability template is fixed for symbols between 1 and 
% 5000.

%--------------------------- ARITHMETIC DECODER ---------------------------
% 'lossless_data': binary arithmetic code
% 'seq': sequence of symbols to be decoded

dseq = zeros(frame_length, n_frames);
dflat_data = zeros(frame_length, n_frames);
    
%------------------------- PROBABILITY TEMPLATE -----------------------
% The probability template is a table containing a set of probability 
% density values, which are trained from a large amount of audio data.
s = 1:1:5000;
prob_template = round(1e4*gaussmf(s,[0.3*max(max(s)) -0.1]));

for i = 1:n_frames
    % Decoding
    dseq(:,i) = arithdeco(lossless_data{1,i}, prob_template, frame_length);

    % Undo avoid zero values by adding 1
    dflat_data(:,i) = dseq(:,i) - 1 ;
    
end

% isequal(seq,dseq) % Check that dseq matches the original seq.

end

