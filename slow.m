function separateUsingSVD(sampleNumber)
    sampleNumber = strcat(num2str(sampleNumber), '/');

    mixedSourcesDir = strcat('sound_files/mixed_sources/', sampleNumber);
    originalSourcesDir = strcat('sound_files/original_sources/', sampleNumber);
    outputDir = strcat('sound_files/svd_output/', sampleNumber);

    numSamples = size(readdir(mixedSourcesDir), 1) - 2;

    [mixedSignal, fs] = audioread(strcat(mixedSourcesDir, '1.wav'));
    [originalSignal, fs] = audioread(strcat(originalSourcesDir, '1.wav'));

    mixedSignalList = [mixedSignal];
    originalSignalList = [originalSignal];

    for i = 2:numSamples
        [mixedSignal, fs] = audioread(strcat(mixedSourcesDir, num2str(i), '.wav'));
        [originalSignal, fs] = audioread(strcat(originalSourcesDir, num2str(i), '.wav'));

        mixedSignalList(:, i) = mixedSignal;
        originalSignalList(:, i) = originalSignal;
    endfor

    xx = mixedSignalList';
    yy = sqrtm(inv(cov(xx')))*(xx-repmat(mean(xx,2),1,size(xx,2)));
    [W,s,v] = svd((repmat(sum(yy.*yy,1),size(yy,1),1).*yy)*yy');

    a = W*xx;

    numGraph = 1;
    for i = 1:numSamples
        subplot(3, numSamples, numGraph);
        plot(originalSignalList(:, i));
        title(strcat('Original wave', ' ', num2str(i)));

        subplot(3, numSamples, numGraph + numSamples);
        plot(mixedSignalList(:, i));
        title(strcat('Mixed wave', ' ', num2str(i)));

        subplot(3, numSamples, numGraph + numSamples * 2);
        plot(a(i, :));
        title(strcat('Output wave', ' ', num2str(i)));

        audiowrite(strcat(outputDir, num2str(i), '.wav'), a(i, :), fs);

        numGraph += 1;
    endfor
endfunction

separateUsingSVD(1);

pause()
