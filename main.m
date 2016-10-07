function [ ] = main( audio_signal )
% MAIN Lossless audio compression tool, which utilizes a pre-processing 
% procedure for flattening the amplitude envelop of linear prediction 
% residue, and an arithmetic coder that adopts a scaled probability template. 

%Play audio input signal
[audio_input,Fs] = audioread(audio_signal);
soundsc(audio_input, Fs);
pause(length(audio_input)/Fs);

% ------------------------Predictor step-----------------------------------
[residues, signs, q_parq] = predictor( audio_signal );

% -----------------------PRE-PROCESSOR STEP--------------------------------
[flat_data, LSB_data] = pre_processor( q_parq, residues );

% -----------------------ENTROPY-ENCODER STEP------------------------------
lossless_data  = entropy_encoder( flat_data );

% -----------------------ENTROPY-DECODER STEP------------------------------
dflat_data = entropy_decoder( lossless_data );

% ------------------------POST-PROCESSOR STEP------------------------------
post_residues = post_processor( signs, dflat_data, LSB_data, q_parq );

% -------------------------RECONSTRUCTOR STEP------------------------------
audio_output = reconstructor( q_parq, post_residues);

% Play and write audio output signal
soundsc(audio_output, Fs);
audiowrite('audio_output.wav',audio_output,Fs);

end

