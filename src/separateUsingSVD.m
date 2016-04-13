pkg load nan
graphics_toolkit gnuplot

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

    st = mktime(localtime(time()));
    xx = mixedSignalList';
    yy = sqrtm(inv(cov(xx')))*(xx-repmat(mean(xx,2),1,size(xx,2)));
    [W,s,v] = svd((repmat(sum(yy.*yy,1),size(yy,1),1).*yy)*yy');
    et = mktime(localtime(time()));

    'TIME'
    extime = et - st
    extime

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

    % f = figure('visible','off');
    % scatter(originalSignalList(:, 1), originalSignalList(:, 2), 1);
    % title('Original signal distribution');
    % xlabel('Original Signal 1')
    % ylabel('Original Signal 2')
    % filename=sprintf(strcat(outputDir,'origDist.png'),f);
    % print(filename);

    % f = figure('visible','off');
    % scatter(mixedSignalList(:, 1), mixedSignalList(:, 2), 1);
    % title('mixed signal distribution');
    % xlabel('Mixed Signal 1')
    % ylabel('Mixed Signal 2')
    % filename=sprintf(strcat(outputDir,'mixDist.png'),f);
    % print(filename);

    % f = figure('visible','off');
    % scatter(outputList(:, 1), outputList(:, 2), 1);
    % title('Output signal distribution');
    % xlabel('Output Signal 1')
    % ylabel('Output Signal 2')
    % filename=sprintf(strcat(outputDir,'outDist.png'),f);
    % print(filename);

    c1 = corrcoef(abs(originalSignal1), abs(output1));
    c2 = corrcoef(abs(originalSignal2), abs(output2));
    cf = mean([c1; c2]);

    c1 = corrcoef(abs(originalSignal1), abs(output2));
    c2 = corrcoef(abs(originalSignal2), abs(output1));
    cs = mean([c1; c2]);

    'CORRELATION:'
    cor = max([cf; cs])
    cor

    jsonString = strcat('{"time":', num2str(extime), ',"cor":', num2str(cor), '}');
    strcat(outputDir, 'stats.json')
    fid = fopen(strcat(outputDir, 'stats.json'), "w")
    fputs(fid, jsonString)

endfunction

function singleSineSVD()
    'INSIDE SINE SVD'
    sample = 'singleSine'
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
    extime = et - st
    extime

    a = W*xx;

    output1 = a(1, :)';
    output2 = a(2, :)';

    mixedSignalList = S;
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

    % f = figure('visible','off');
    % scatter(originalSignalList(:, 1), originalSignalList(:, 2), 1);
    % title('Original signal distribution');
    % xlabel('Original Signal 1')
    % ylabel('Original Signal 2')
    % filename=sprintf(strcat(outputDir,'origDist.png'),f);
    % print(filename);

    % f = figure('visible','off');
    % scatter(mixedSignalList(:, 1), mixedSignalList(:, 2), 1);
    % title('mixed signal distribution');
    % xlabel('Mixed Signal 1')
    % ylabel('Mixed Signal 2')
    % filename=sprintf(strcat(outputDir,'mixDist.png'),f);
    % print(filename);

    % f = figure('visible','off');
    % scatter(outputList(:, 1), outputList(:, 2), 1);
    % title('Output signal distribution');
    % xlabel('Output Signal 1')
    % ylabel('Output Signal 2')
    % filename=sprintf(strcat(outputDir,'outDist.png'),f);
    % print(filename);

    c1 = corrcoef(abs(originalSignal1), abs(output1));
    c2 = corrcoef(abs(originalSignal2), abs(output2));
    cf = mean([c1; c2]);

    c1 = corrcoef(abs(originalSignal1), abs(output2));
    c2 = corrcoef(abs(originalSignal2), abs(output1));
    cs = mean([c1; c2]);

    'CORRELATION'
    cor = max([cf, cs])
    cor

    jsonString = strcat('{"time":', num2str(extime), ',"cor":', num2str(cor), '}');
    fid = fopen(strcat(outputDir, 'stats.json'), "w")
    fputs(fid, jsonString)
endfunction

sample = argv(){1,1}

if (strcmp(sample, 'singleSine'))
    singleSineSVD()
else
    separateUsingSVD(sample)
endif

pause()
