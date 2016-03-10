function separateUsingSVD(sampleNumber)
    sampleNumber = strcat(num2str(sampleNumber), '/');

    mixedSourcesDir = strcat('sound_files/mixed_sources/', sampleNumber);
    originalSourcesDir = strcat('sound_files/original_sources/', sampleNumber);
    outputDir = strcat('sound_files/svd_output/', sampleNumber);

    numSamples = size(readdir(mixedSourcesDir), 1) - 2;

    [mixedSignal1, fs] = audioread(strcat(mixedSourcesDir, '1.wav'));
    [originalSignal1, fs] = audioread(strcat(originalSourcesDir, '1.wav'));
    [mixedSignal2, fs] = audioread(strcat(mixedSourcesDir, '2.wav'));
    [originalSignal2, fs] = audioread(strcat(originalSourcesDir, '2.wav'));

    mixedSignalList = [mixedSignal1, mixedSignal2];
    originalSignalList = [originalSignal1, originalSignal2];

    xx = mixedSignalList';
    yy = sqrtm(inv(cov(xx')))*(xx-repmat(mean(xx,2),1,size(xx,2)));
    st = mktime(localtime(time()));
    [W,s,v] = svd((repmat(sum(yy.*yy,1),size(yy,1),1).*yy)*yy');
    et = mktime(localtime(time()));

    'TIME'
    et - st

    a = W*xx;

    output1 = a(1, :)';
    output2 = a(2, :)';
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

        audiowrite(strcat(outputDir, num2str(i), '.wav'), a(i, :), fs);

        numGraph += 1;
    endfor

    originalSignal1 = 255 * originalSignal1 / norm(originalSignal1);
    originalSignal2 = 255 * originalSignal2 / norm(originalSignal2);
    output1 = 255 * output1 / norm(output1);
    output2 = 255 * output2 / norm(output2);

    e1 = abs((abs(originalSignal1) - abs(output1)));
    e1(! isfinite(e1)) = 0;
    e1 = sqrt(mean(e1 .^ 2));
    e2 = abs((abs(originalSignal2) - abs(output2)));
    e2(! isfinite(e2)) = 0;
    e2 = sqrt(mean(e2 .^ 2));
    ef = abs(mean([e1 ./ mean(abs(originalSignal1)); e2 ./ mean(abs(originalSignal2))]));

    e1 = abs((abs(originalSignal1) - abs(output2)));
    e1(! isfinite(e1)) = 0;
    e1 = sqrt(mean(e1 .^ 2));
    e2 = abs((abs(originalSignal2) - abs(output1)));
    e2(! isfinite(e2)) = 0;
    e2 = sqrt(mean(e2 .^ 2));
    es = abs(mean([e1 ./ mean(abs(originalSignal1)); e2 ./ mean(abs(originalSignal2))]));
    'ERROR'
    min([ef, es])
endfunction

function customMixAndSeparate(sample)
    sample = strcat(num2str(sample), '/');

    mixedSourcesDir = strcat('sound_files/mixed_sources/', sample);
    originalSourcesDir = strcat('sound_files/original_sources/', sample);
    outputDir = strcat('sound_files/svd_output/', sample);

    [originalSignal1, fs] = audioread(strcat(originalSourcesDir, '1.wav'));
    [originalSignal2, fs] = audioread(strcat(originalSourcesDir, '2.wav'));

    A = unifrnd(1, 5, 2, 2);
    M = [originalSignal1, originalSignal2];
    S = M * A;
    xx = S';
    yy = sqrtm(inv(cov(xx')))*(xx-repmat(mean(xx,2),1,size(xx,2)));
    st = mktime(localtime(time()));
    [W,s,v] = svd((repmat(sum(yy.*yy,1),size(yy,1),1).*yy)*yy');
    et = mktime(localtime(time()));

    'TIME'
    et - st

    a = W*xx;

    output1 = a(1, :)';
    output2 = a(2, :)';

    mixedSignalList = M;
    outputList = [output1, output2];
    originalSignalList = [originalSignal1, originalSignal2];

    numGraph = 1;
    numSamples = 2;
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

        audiowrite(strcat(outputDir, num2str(i), '.wav'), a(i, :), fs);

        numGraph += 1;
    endfor

    originalSignal1 = 255 * originalSignal1 / norm(originalSignal1);
    originalSignal2 = 255 * originalSignal2 / norm(originalSignal2);
    output1 = 255 * output1 / norm(output1);
    output2 = 255 * output2 / norm(output2);

    e1 = abs((abs(originalSignal1) - abs(output1)));
    e1(! isfinite(e1)) = 0;
    e1 = sqrt(mean(e1 .^ 2));
    e2 = abs((abs(originalSignal2) - abs(output2)));
    e2(! isfinite(e2)) = 0;
    e2 = sqrt(mean(e2 .^ 2));
    ef = abs(mean([e1 ./ mean(abs(originalSignal1)); e2 ./ mean(abs(originalSignal2))]));

    e1 = abs((abs(originalSignal1) - abs(output2)));
    e1(! isfinite(e1)) = 0;
    e1 = sqrt(mean(e1 .^ 2));
    e2 = abs((abs(originalSignal2) - abs(output1)));
    e2(! isfinite(e2)) = 0;
    e2 = sqrt(mean(e2 .^ 2));
    es = abs(mean([e1 ./ mean(abs(originalSignal1)); e2 ./ mean(abs(originalSignal2))]));
    'ERROR'
    min([ef, es])
endfunction

customMixAndSeparate('doubleSine');

% separateUsingSVD(5)

pause()
