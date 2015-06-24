        curSegment = (i-1)*step+1 : (i-1)*step+window;
        prominence = 85;
        RAW = PPG(1,index).sig;
        % MSSA Using both ppg channels.
        RAW1 = RAW(2,:);
        RAW2 = RAW(3,:);
        BPM0 = BPM(1,index).BPM0;
        display(numel(BPM0));
        display(index);
        modPPG1 = RAW1;
        modPPG2 = RAW2;
        for point = 1:numel(RAW1)
            if RAW1(point) <= 0
                modPPG1(point) = 0;
            end
        end
        for point = 1:numel(RAW2)
            if RAW2(point) <= 0
                modPPG2(point) = 0;
            end
        end
        clearvars point;
        [peaks1,locations1] = findpeaks(RAW1);
        [peaksFilt1,locationsFilt1] = findpeaks(RAW1,'MinPeakDistance',prominence);
        [peaks2,locations2] = findpeaks(RAW2);
        [peaksFilt2,locationsFilt2] = findpeaks(RAW1,'MinPeakDistance',prominence);

        locationCount1 = numel(locations1);
        locationCount2 = numel(locations2);
        for count = 1:locationCount1
            element = locations1(count);
            if count < locationCount1
                nextElement = locations1(count+1);
            else
                nextElement = element;
            end
            nextDifferenceSeconds1(count) = (nextElement - element)/125; %div 125 converts to seconds
        end

        for count = 1:locationCount2
            element = locations2(count);
            if count < locationCount2
                nextElement = locations2(count+1);
            else
                nextElement = element;
            end
            nextDifferenceSeconds2(count) = (nextElement - element)/125; %div 125 converts to seconds
        end
        
        count1 = 1;
        n1 = 1;
        maximum1 = numel(RAW1);
        maximum2 = numel(locationsFilt1);
        while count1 < maximum1
            if n1 < maximum2
                if count1 == locationsFilt1(n1)
                    modPPG1(count1) = peaksFilt1(n1);
                    n1 = n1 + 1;
                end
            end
            count1 = count1 + 1;
        end
        start1 = numel(modPPG1);
        maximumA = numel(RAW1);
        if start1 < maximumA
            for count1 = start1:maximumA
                modPPG1(count1) = RAW1(count1);
            end
        end
        resultSawtelle1 = modPPG1;
        
        count2 = 1;
        n2 = 1;
        maximum1 = numel(RAW2);
        maximum2 = numel(locationsFilt2);
        while count2 < maximum1
            if n2 < maximum2
                if count2 == locationsFilt2(n2)
                    modPPG2(count2) = peaksFilt2(n2);
                    n2 = n2 + 1;
                end
            end
            count2 = count2 + 1;
        end
        start2 = numel(modPPG2);
        maximumB = numel(RAW2);
        if start2 < maximumB
            for count2 = start2:maximumB
                modPPG2(count2) = RAW2(count2);
            end
        end
        resultSawtelle2 = modPPG2;
        windowSize=1000;
        stepSize=625; % robut for 250-625, but gets noticeably slower with smaller
                      % step size.
        L=floor(windowSize/2);
        K=windowSize-L+1;
        N=4096;
        fs = 125;
        originalLength1 = length(modPPG1);
        originalLength2 = length(modPPG2);
        winNum1 = floor(originalLength1/stepSize);
        winNum2 = floor(originalLength2/stepSize);
        winNum = min(winNum1,winNum2);
        tempPpg1 = zeros(1,(winNum)*stepSize);
        tempPpg1(1:originalLength1) = modPPG1(1:originalLength1);
        tempPpg2 = zeros(1,(winNum)*stepSize);
        tempPpg2(1:originalLength2) = modPPG2(1:originalLength2);
        windowedSignal1 = zeros(winNum, windowSize);
        windowedSignal2 = zeros(winNum, windowSize);
        for i = 0:winNum-1 % 0-index for first row with no delay (aka window starting point)
            for j = 1:windowSize
                if (j+(stepSize*i) > originalLength1) 
                    windowedSignal1(i+1, j) = 0;
                    windowedSignal2(i+1, j) = 0;
                else
                    windowedSignal1(i+1,j) = tempPpg1(1, j + stepSize*i);
                    windowedSignal2(i+1,j) = tempPpg2(1, j + stepSize*i);
                end
            end
        end
        restoredWindows = zeros(winNum,windowSize);
        for i = 1:winNum
            currentWindow = groupingMSSA(samples1(i,:), samples2(i,:));
            currentWindow = sgolayfilt(currentWindow, 0 , 19);
            currentWindow = currentWindow-mean(currentWindow);
            currentWindow = currentWindow/max(currentWindow);
            restoredWindows(i,:) = currentWindow;
        end
        resultCobo = zeros(1,winNum*stepSize);
        for i = 0:winNum-1
            currentWindow = restoredWindows(i+1, :);
            resultCobo(1 + i*stepSize : windowSize + i*stepSize ) = currentWindow;
        end
        if length(resultCobo) > start1
            resultCobo = resultCobo(1:start1);
        end
        error_og = zeros;
        error_ac = zeros;
        fs=125;
        fmin=1;
        fmax=3.3;
        x1=resultCobo;
        peaks =zeros;
        yg=zeros;
        peaks_ac =zeros;
        collp = zeros;
        error_peaks = zeros;
        peaks1 = zeros; 
        track =zeros; 
        error_peaks1 = zeros;
        max_aa = zeros; 
        L=length(x1(1,:))-(8*fs); %Because first average point is after 8 seconds.
        N=(L/(2*fs)); %number of 2 second periods left. 
        firstrun=1;
        datasetfailed = 0;

        while N > 1 %Need it to be 1 to go to end without error?
            if firstrun ==1; 
                x2=x1(1,1:(8*fs));
                last=8*fs;
                firstrun = 0;
                n=1;
            else 
                last=last+2*fs;
                x2=x1(1,(last-8*fs):last);

                %N = N-2; wtf??? This was a big mistake.... ( .__.)
                N=N-1;
                n=n+1;
            end
            ac=xcorr(x2);
            t1=0:(1/125):(length(ac)-1)/125;
            t=0:(1/125):(length(x2)-1)/125;
            NFFT=4096;
            f=fs/2*linspace(0,1,NFFT/2+1);
            Y=fft(ac, NFFT)/(length(x2)); % Take fft take
            Y=2*abs(Y(1:NFFT/2+1)); %Take care of 2 factor, and absolute
            roii=fmax>f & f>fmin; %Setup frequency range that is realistic as well. 
            Y=Y(roii)/max(Y(roii)); %Set specific region up, and normalize fft. 
            f1=f(roii);
            roii=fmax>f & f>fmin;
            delta=6;
            pdm = 1; %default: peak loctioan minus delta equals zero. 
            pdp = length(f1); %default: peak loctioan plus delta equals end of spectrum
            win = 0.05;
            if (n>4)
                if (max_aa(end) < win & max_aa(end-1) < win)
                    delta=delta+2;
                end
                if (max_aa(end) < win & max_aa(end-1) < win & max_aa(end-2) < win)
                    delta=delta+4;
                end
            end
            if (n<2)
                roi=1.98>f1 & f1>1.0;
            else
                pdp=windowMaxFFT+delta;
                pdm=windowMaxFFT-delta;
            end

            if(pdp > length(f1))
                pdm=pdm-(abs(pdp)-length(f1));
                pdp = length(f1);
            end
            if (pdm <1)
                pdp=pdp+abs(pdm);
                if (pdp >length(f1))
                    pdp = length(f1);
                end
                pdm =1;
            end
            if (n>2)
                pdp;
                pdm;
                roi = f1(pdp)>f1 & f1>f1(pdm);
            end
            windowMaxFFT=find(Y(:) == max(Y(roi)));
            max_a=Y(windowMaxFFT);
            max_aa=[max_aa max_a];
            actual=BPM0(n,1);
            track(n)=1;
            peaks1 = [peaks1, 60*f1(windowMaxFFT)];
            flag20= abs(peaks1(end)-actual)
            if (flag20 > 20)
                datasetfailed = 1;
            end
            peaksize=length(peaks1);
            error_acp=(abs(60*f1(windowMaxFFT)-BPM0(n,1))/BPM0(n,1))*100;
            error_ac=[error_ac error_acp];
            peaks_ac= [peaks_ac abs(60*f1(windowMaxFFT))];
            error_peaks =[error_peaks ((peaks(end)-BPM0(n,1))/BPM0(n,1))*100];
            error_peaks1 =[error_peaks1 abs(((peaks1(end)-BPM0(n,1))/BPM0(n,1))*100)];
            hold off;pause(.0001);
        end
        peaks1=smooth(peaks1,0.05, 'rloess');
        avg_peaks1_error = (sum(error_peaks1)/length(error_peaks1));
        avg_error_ac=(sum(error_ac)/length(error_ac));
        Error(index,2) = avg_error_ac;
