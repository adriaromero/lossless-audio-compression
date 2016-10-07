function [residues, signs, q_parq] = predictor( audio_input )
% PREDICTOR Input audio samples are segmented into frames of fixed 
% length. Linear predictive coding (LPC) is then performed on each frame, 
% with partial-correlation (PARCOR) coefficients computed through the 
% Levinson-Durbin algorithm. The PARCOR coefficients are quantized and sent
% as ancillary information in the lossless bitstream. 
% The quantized PARCOR coefficients are also locally dequantized and con-
% verted to the tap coefficients of a linear predictor, which generates a 
% prediction for each sample in the frame. 
% The difference signal between an input sample and its prediction - 
% the prediction residue - is output to the next processing stage.
%
% INPUT VARIABLES: 
%       'audio_input': Audio input signal with format .wav, .ogg, .flac,
%       .mp3 or .mp4
%
% OUTPUT VARIABLES:
%       'Fs': Sample rate, in Hz, of audio_input.
%       'residues': Prediction residues
%       'signs': Sign (-1,0,1) of each prediction residues values
%       'q_parq':   Quantized PARCOR coefficients. 
%                   Matrix dimensions: (parcor_order,n_frames)

%----------------------------READING INPUT FILE----------------------------
% 'audio': Audio data in the file, returned as an m-by-n matrix, where m is 
% the number of audio samples read and n is the number of audio channels in 
% the file.
% 'Fs': Sample rate, in Hz, of audio data 'y', returned as a positive scalar.
%'n_samples': Number of total samples of audio data
audio_input = 'audio_input.wav';
[audio,Fs] = audioread(audio_input);
audio = audio * 2.^15;  % original values into integer
n_samples = length(audio);

%----------------------------FRAMING---------------------------------------
% Objective: to segment the audio_input data into frames of fixed length 
% 'frame_length': fixed length of the audio frame
% 'n_frames': Number of frames
% 'frames': Array that contains all frames data
%           Dimension = (frame_length, n_frames)
%           i.e. frames(67,3) : sample 67 of frame 3

frame_length = Fs/40;    % Assuming Fs is integer. 
                         % Frame length of 400 ms (for Fs = 16kHz) 
trailing_samples = mod(n_samples, frame_length);
frames = reshape( audio(1:end-trailing_samples), frame_length, []);
n_frames = length(frames(1,:));

    
%----------------------------PARCOR----------------------------------------
% ACF	Autocorrelation function from lag=[0:p]
% [a,e,k] = levinson(r,n) returns the coefficients for an autoregressive 
% model of order n (a), the prediction error (e) of order n and the 
% reflection coefficients k as a column vector of length n.
% 'parq': PArtian AutoCORrelation Coefficients 
% 'RC': Reflection Coefficient

parcor_order = 20; 
parq = zeros(parcor_order,n_frames);
lpc_original = zeros(parcor_order+1,n_frames);
wind = hamming(frame_length);   %Hamming windowing

for i = 1:n_frames
    ACF = autocorr(frames(:,i) .* wind);
    [a,e,k] = levinson(ACF,parcor_order);
    lpc_original(:,i) = a;
    parq(:,i)= k;    
end

%----------------------------PARCOR QUANTIZATION---------------------------
% Parcor coefficients (parq, k = 1,...,lpc_order) are first quantized by 
% the companding function. The resulting quantized values, 'q_parq', are 
% restricted to the range [-64,63].

q_parq = zeros(parcor_order,n_frames);
for i = 1:n_frames
    q_parq(1,i) = floor(64*log(2/3 + 5/6*sqrt((1+parq(1,i))/2))/log(3/2)); 
    q_parq(2,i) = floor(64*log(2/3 + 5/6*sqrt((1-parq(2,i))/2))/log(3/2));
    q_parq(3:parcor_order,i) = floor(64*parq(3:parcor_order,i));
end

%----------------------------PARCOR DE-QUANTIZATION------------------------
% 'deq_parq': De-quantizated PARCOR coefficients
deq_parq = zeros(parcor_order,n_frames);

for i = 1:n_frames
    deq_parq(1,i) = (2 * ((exp(q_parq(1,i)/64*log(3/2))-(2/3)) * 6/5).^2)-1; 
    deq_parq(2,i) = -(2 * ((exp(q_parq(2,i)/64*log(3/2))-(2/3)) * 6/5).^2)+1;
    deq_parq(3:parcor_order,i) = q_parq(3:parcor_order,i)/64; 
end

%----------------------------PARCOR TO LPC---------------------------------
% The reconstructed Parcor coefficients are converted to k-order 
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

%----------------------------LINEAR PREDICTOR ENCODER----------------------
% The Linear Predictor generates a prediction for each sample in the frame.
% After that, it computes the prediction residues.
% 'residues': Prediction residues

residues = zeros(parcor_order,n_frames);

for i = 1:n_frames
    residues(1,i) = frames(1,i);
    
    for j = 2:parcor_order
        sum = 0;
        for k = 1:j-1
            sum = sum + lpc(k,j,i)*frames((j-k),i);
        end 
        residues(j,i) = frames(j,i)+round(sum);
    end
   
    for j = parcor_order+1:frame_length
        sum = 0;
        for k = 1:parcor_order
            sum = sum + lpc(k,parcor_order,i)*frames(j-k,i);
        end 
        residues(j,i) = frames(j,i)+round(sum);
    end
end

%---------------------------SIGN VALUE-------------------------------------
signs = sign(residues);

end

