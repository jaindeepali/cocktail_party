clear;
[x1, f] = wavread('sound_files/mixed_sources/4/1.wav');
x2 = wavread('sound_files/mixed_sources/4/2.wav');
x = [x1 x2];
[rx,cx] = size(x);
T = 512;
Q = 128;
K = 5;
N = floor(rx/T/K);
verbose = 0;
iter = 1000;
printevery = 10;
eta = 1.0;
ds =2;

[y,Wt,J] = convbss(x,T,Q,K,N,verbose,iter,printevery,eta,ds);

y1 = y(:,1);
y2 = y(:,2);
wavwrite( y1, f, 16, 'sound_files/conv_output/4/1.wav' );
wavwrite( y2, f, 16, 'sound_files/conv_output/4/2.wav' );