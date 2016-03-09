[orig1, _1] = audioread('sound files/original sources/source1.wav');
[orig2, _2] = audioread('sound files/original sources/source2.wav');

[x1, Fs1] = audioread('sound files/mixed sources/X1_linear.wav');
[x2, Fs2] = audioread('sound files/mixed sources/X2_linear.wav');
xx = [x1, x2]';
yy = sqrtm(inv(cov(xx')))*(xx-repmat(mean(xx,2),1,size(xx,2)));
[W,s,v] = svd((repmat(sum(yy.*yy,1),size(yy,1),1).*yy)*yy');

a = W*xx; %W is unmixing matrix
subplot(3,2,1); plot(x1); title('mixed audio - mic 1');
subplot(3,2,2); plot(x2); title('mixed audio - mic 2');
subplot(3,2,3); plot(a(1,:), 'g'); title('unmixed wave 1');
subplot(3,2,4); plot(a(2,:), 'r'); title('unmixed wave 2');
subplot(3,2,5); plot(orig1); title('original wave 1');
subplot(3,2,6); plot(orig2); title('original wave 2');

audiowrite('unmixed1.wav', a(1,:), Fs1);
audiowrite('unmixed2.wav', a(2,:), Fs1);

pause()
