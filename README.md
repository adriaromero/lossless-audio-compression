# LOSSLESS AUDIO COMPRESSION IN IEEE 1857.2 STANDARD FOR ADVANCED AUDIO CODING

# Authors: Adrià Romero López, Dan-Mihai Bădescu <br>
Universitat Politècnica de Catalunya (UPC), Barcelona, 2015 

# ABSTRACT
In August 2013, IEEE approved a new standard that comprises both lossy and lossless audio compression tools. This new standard is called IEEE 1857.2. This paper focus on the lossless audio compression tool, which utilizes a pre- processing procedure for flattening the amplitude envelop of the linear prediction residue, and an arithmetic encoder that adopts a scaled probability template.
Index Terms—Lossless audio compression for advanced audio coding, IEEE 1857.2

# INTRODUCTION
Generally, the multimedia contents are pre-compressed before they are uploaded online. Normally, there are two types of compression methods: the lossy and the lossless. The lossy method attempts to remove perceptually less important information from the audio data while keeping the sound quality very close to the original one. Some examples that include this type of compression method are: MPEG-1 Layer (MP3) and the MPEG-2/4 Advanced Audio Coding (AAC) which achieve more than twenty times compression and still delivering good sound quality. On the other hand, the lossless method keeps every bit of information from the original audio data and can achieve about two times compression.
Lossy audio compression is mainly used for playing music from iPods, or listening to networked radio using mobile phones. In contrast, lossless audio reproduction, archival of database and more recently biomedical signal compression, such as lossless ECG compression [1].
A generally approach used in lossless audio compression, is a combination of linear prediction and entropy encoding. The linear predictor first removes the redundancy in the input data and generates a prediction residue, which is encoded by the entropy encoder. The system is based on a LPC predictor, a pre-processor for flattening the prediction residue, and an entropy coder that is based on arithmetic coding with probability template scaling [1].

