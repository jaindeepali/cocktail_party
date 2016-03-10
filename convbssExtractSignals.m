function convbssExtractSignals(sampleNumber)
    sampleNumber = strcat(num2str(sampleNumber), '/');

    mixedSourcesDir = strcat('sound_files/mixed_sources/', sampleNumber);
    originalSourcesDir = strcat('sound_files/original_sources/', sampleNumber);
    outputDir = strcat('sound_files/conv_output/', sampleNumber);

    numSamples = size(readdir(mixedSourcesDir), 1) - 2;

    [mixedSignal1, fs] = audioread(strcat(mixedSourcesDir, '1.wav'));
    [originalSignal1, fs] = audioread(strcat(originalSourcesDir, '1.wav'));
    [mixedSignal2, fs] = audioread(strcat(mixedSourcesDir, '2.wav'));
    [originalSignal2, fs] = audioread(strcat(originalSourcesDir, '2.wav'));

    mixedSignalList = [mixedSignal1, mixedSignal2];
    originalSignalList = [originalSignal1, originalSignal2];

    [rx,cx] = size(mixedSignalList);
    T = 512;
    Q = 128;
    K = 5;
    N = floor(rx/T/K);
    verbose = 0;
    iter = 1000;
    printevery = 10;
    eta = 1.0;
    ds =2;
    
    st = mktime(localtime(time()));

    [y,Wt,J] = convbss(mixedSignalList,T,Q,K,N,verbose,iter,printevery,eta,ds);

    et = mktime(localtime(time()));

    'TIME'
    et - st

    output1 = y(:, 1);
    output2 = y(:, 2);
    outputList = [output1, output2];

    numGraph = 1;
    for i = 1:numSamples
        subplot(3, numSamples, numGraph);
        plot(originalSignalList(:, i));
        title(strcat('Original wave', ' ', num2str(i)));

        subplot(3, numSamples, numGraph + numSamples);
        plot(mixedSignalList(:, i));
        title(strcat('Mixed wave', ' ', num2str(i)));

        subplot(3, numSamples, numGraph + numSamples * 2);
        plot(outputList(:, i));
        title(strcat('Output wave', ' ', num2str(i)));

        audiowrite(strcat(outputDir, num2str(i), '.wav'), y(:, 1)', fs);

        numGraph += 1;
    endfor

    c1 = corrcoef(abs(originalSignal1), abs(output1));
    c2 = corrcoef(abs(originalSignal2), abs(output2));
    cf = mean([c1; c2]);

    c1 = corrcoef(abs(originalSignal1), abs(output2));
    c2 = corrcoef(abs(originalSignal2), abs(output1));
    cs = mean([c1; c2]);

    'CORRELATION:'
    max([cf; cs])
endfunction

pause()

convbssExtractSignals(5)