function [ calculatedBPM, windowMaxFFTindexArray, windowMaxFFTAmplitudeArray ] = MACVSSA( PPG1, PPG2,  windowMaxFFTindexArray, windowMaxFFTAmplitudeArray, n )
        prominence = 1;
        RAW1 = PPG1;
        RAW2 = PPG2;
        modPPG1 = PPG1;
        modPPG2 = PPG2;
        roi = 1;
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
        [peaksFilt1,locationsFilt1] = findpeaks(modPPG1,'MinPeakDistance',prominence);
        [peaksFilt2,locationsFilt2] = findpeaks(modPPG2,'MinPeakDistance',prominence);
        
        
        for i = numel(locationsFilt1)
           modPPG1(locationsFilt1) = peaksFilt1; 
        end
        for i = numel(locationsFilt2)
           modPPG2(locationsFilt2) = peaksFilt2; 
        end
        
%         count1 = 1;
%         n1 = 1;
%         maximum1 = numel(RAW1);
%         maximum2 = numel(locationsFilt1);
%         while count1 < maximum1
%             if n1 < maximum2
%                 if count1 == locationsFilt1(n1)
%                     modPPG1(count1) = peaksFilt1(n1);
%                     n1 = n1 + 1;
%                 end
%             end
%             count1 = count1 + 1;
%         end
%         start1 = numel(modPPG1);
%         maximumA = numel(RAW1);
%         if start1 < maximumA
%             for count1 = start1:maximumA
%                 modPPG1(count1) = RAW1(count1);
%             end
%         end
%         
%         count2 = 1;
%         n2 = 1;
%         maximum1 = numel(RAW2);
%         maximum2 = numel(locationsFilt2);
%         while count2 < maximum1
%             if n2 < maximum2
%                 if count2 == locationsFilt2(n2)
%                     modPPG2(count2) = peaksFilt2(n2);
%                     n2 = n2 + 1;
%                 end
%             end
%             count2 = count2 + 1;
%         end
%         start2 = numel(modPPG2);
%         maximumB = numel(RAW2);
%         if start2 < maximumB
%             for count2 = start2:maximumB
%                 modPPG2(count2) = RAW2(count2);
%             end
%         end
%         figure
%         plot(modPPG1);
%         figure
%         plot(modPPG2);

%         % Was causing length error when error.
%         modPPG1 = sgolayfilt(modPPG1(1,:), 0, 19);
%         modPPG2 = sgolayfilt(modPPG2(1,:), 0, 19);

       MPPG1=modPPG1(1,:);
       MPPG2=modPPG2(1,:);
%        plot(MPPG1)
%        pause(0.1);
         % Testing introducing bp filt again.
         fs=125;
        OR=6; %Order
        fc1= 1; %1st Frequency cutoff
        fc2 =3.33;%2nd Freq. cutoff
        [b,a] = butter(OR, [fc1 fc2]/(fs/2), 'bandpass');
        MPPG1 = filter(b,a,MPPG1);
        MPPG2 = filter(b,a,MPPG2);
%        plot(MPPG1)
%        pause(0.1);
        currentWindow = groupingMSSA(MPPG1, MPPG2);
        length(currentWindow)
        currentWindow = sgolayfilt(currentWindow, 0 , 19);
        length(currentWindow)

        currentWindow = currentWindow-mean(currentWindow);
        length(currentWindow)
        currentWindow = currentWindow/max(currentWindow);
        length(currentWindow)
        fs=125;
        fmin=1;
        fmax=3.3;
        
        % Only need current window.
        x=currentWindow; 
        
        % Testing introducing bp filt again.
        OR=6; %Order
        fc1= 1.4; %1st Frequency cutoff
        fc2 =2.3;%2nd Freq. cutoff
        [b,a] = butter(OR, [fc1 fc2]/(fs/2), 'bandpass');
        x1 = filter(b,a,x);
        
        x2 = x1;
        
        ac=xcorr(x2);
        NFFT=4096;
        f=fs/2*linspace(0,1,NFFT/2+1);
        Y=fft(ac, NFFT)/(length(x2)); % Take fft take
        Y=2*abs(Y(1:NFFT/2+1)); %Take care of 2 factor, and absolute
        roii=fmax>f & f>fmin; %Setup frequency range that is realistic as well. 
        Y=Y(roii)/max(Y(roii)); %Set specific region up, and normalize fft. 
        f1=f(roii);
        delta=6;
        pdm = 1; %default: peak loctioan minus delta equals zero. 
        pdp = length(f1); %default: peak loctioan plus delta equals end of spectrum
        win = 0.05;
        
        if (n>4)
            if (windowMaxFFTAmplitudeArray(end) < win && windowMaxFFTAmplitudeArray(end-1) < win)
                delta=delta+2;
            end
            if (windowMaxFFTAmplitudeArray(end) < win && windowMaxFFTAmplitudeArray(end-1) < win && windowMaxFFTAmplitudeArray(end-2) < win)
                delta=delta+4;
            end
        end
        if (n<2)
            roi=1.98>f1 & f1>1.0;
        else
            pdp=windowMaxFFTindexArray(end)+delta;
            pdm=windowMaxFFTindexArray(end)-delta;
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
            roi = f1(pdp)>f1 & f1>f1(pdm);
        end
        windowMaxFFTindex=find(Y(:) == max(Y(roi)));
        windowMaxFFTAmplitude=Y(windowMaxFFTindex);
        windowMaxFFTAmplitudeArray = [windowMaxFFTAmplitudeArray, windowMaxFFTAmplitude];        
        windowMaxFFTindexArray = [windowMaxFFTindexArray, windowMaxFFTindex];
        calculatedBPM=60*f1(windowMaxFFTindex);
        %display(calculatedBPM);
        
        %Trying to visualizing something.
        subplot(1,2,1),plot(ac);
        subplot(1,2,2),plot(Y);
        hold on;
        subplot(1,2,2),stem(windowMaxFFTindex,windowMaxFFTAmplitude);
        hold off;
        pause(0.01)
end

