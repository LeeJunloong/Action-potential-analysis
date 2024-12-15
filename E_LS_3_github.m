% 20220420, F-I curve, complexity spike (based on half-width) and sodium spike (based on slope) finder
% by LJL
clear all;clc
Fidx{1}=dir('test.mat'); % filename
slopesodium= 20000; % mV/s, can be modified
% slopecalcium= 2000; % mV/s, can be modified. Possibly add an upper limit?
%%
for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        %% loading data
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        pause(1)
        %% ------------- reading data ---------------------------
        clear(('*Ch31*'));
        varia = who('*Ch*')
        for i=1:length(varia)
            if isfield(eval(varia{i}),'units')
                str=eval([varia{i},'.units']);
                if ismember('pA',str)
                    Im=eval([varia{i} '.values']);
                    Im=Im-mean(Im(1:10));
                    dtI=eval([varia{i} '.interval']);
                end
                if ismember('mV',str)
                    Vm=eval([varia{i} '.values']);
                    dtV=eval([varia{i} '.interval']);
                end
            end
        end
        clear (varia{:}) 
       %% find key time points
        TimeStim = currenteventfinder(Im,5,1/dtI,0.2,0.5,1,0).*dtI; % (Im,10,1/dtI,0.2,0.5,1,0) The 4th parameter can control the stimulus duration
        Dura=round(range(TimeStim,2)*100)./100;
        TimeStim=TimeStim(:,1);
        disp(length(TimeStim));
        %% waveform extraction
        pretime=0.1;
        postime=0.8;
        waveformV=[];
        waveformI=[];
        Vrest=[];
        for i=1:length(TimeStim)
            dataidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV+postime/dtV);
            baseidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV);
            dataidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI+postime/dtI);
            baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
            Vrest=[Vrest,mean(Vm(baseidxV))];
            waveformV=[waveformV,Vm(dataidxV)];
            waveformI=[waveformI,Im(dataidxI)-mean(Im(baseidxI))];
        end
        tspanV=linspace(-pretime,postime,size(waveformV,1));
        tspanI=linspace(-pretime,postime,size(waveformI,1));
        wrongidx=[];
        waveformI(:,wrongidx)=[];
        waveformV(:,wrongidx)=[];
        Vrest(wrongidx)=[];
        TimeStim(wrongidx)=[];
        Dura(wrongidx)=[];

        %% create necessary tags
        % spike shape
        spkFreq=[];
        spkTime=[];
        CspkFreq=[];
        CspkTime=[];
        slopwaveformV=[];
        cpeaktime=[];
        sodiumpeaktime=[];
        calciumpeaktime=[];
        WaveformV=[];
        ALLpeakvalue=[];
        ALLpeaktime=[];
        
        for i=1:size(waveformV,2)    
           slopwaveformV = diff(waveformV(:,i))/dtV
           slopwaveformV = smooth(slopwaveformV)
           [Speakvalue,sodiumpeaktime] = findpeaks(slopwaveformV,'minpeakheight',slopesodium,'minpeakdistance',0.003/dtV)
           % Remove tailing stimulus peaks
           sodiumpeaktime(sodiumpeaktime*dtV>(Dura(i)+pretime-0.001)|sodiumpeaktime*dtV<(pretime+0.0002))=[] 
           
           % For each sodium spike, find the sodium action potential peak
           speaktime=[]
           for m=1:length(sodiumpeaktime)
             % find the first peak after threshold
             datapeaks = waveformV(sodiumpeaktime(m,1):(sodiumpeaktime(m,1)+0.005/dtV),i)
             datapeaks=smooth(datapeaks)
             [sopeakvalue,sopeaktime] = findpeaks(datapeaks,'minpeakheight',-10,'minpeakdistance',0.002/dtV)
             if ~isempty(sopeaktime) & datapeaks(sopeaktime)-min(datapeaks)>=20 % amplitude check
                 sopeaktime=sopeaktime
             else
                 sopeaktime=[]
             end
                 
             if ~isempty(sopeaktime)
                sopeaktime=sopeaktime+sodiumpeaktime(m,1)
             end
             speaktime=[speaktime,sopeaktime'] % add new value to the sodium peak set
           end
            spkFreq(i)=length(speaktime)/0.5;
            spkTime{i}=speaktime; 
            WaveformV = smooth(waveformV(:,i))
            % Find slow peaks
            ALLpeaktime = peakfinder(waveformV(:,i)-Vrest(i),10,-10,1,0);
            % AHP = peakfinder(waveformV(:,i)-Vrest(i),10,-20,-1,0)
            % Remove tailing stimulus peaks
            ALLpeaktime(ALLpeaktime*dtV>(Dura(i)+pretime-0.008)| ALLpeaktime*dtV<(pretime+0.0005))=[]
            AHP=[];
            AHPposition=[];
            AHPpos=[];
            AHP1piont=[];
            AHP1pos=[];
            Decaypos=[];
            Risepos=[];
            CpeakHW=[];
            firstHW=[];
            
            if  ~isempty(ALLpeaktime)&&length(ALLpeaktime)==1 % If there is only one peak
                j=length(ALLpeaktime);
                [AHP(j),AHPposition]=min(waveformV(ALLpeaktime(j):(ALLpeaktime(j)+0.05/dtV),i));
                AHPpos(1,j)=AHPposition+ALLpeaktime(j);
                [AHP1piont,AHP1pos]=min(abs((waveformV((ALLpeaktime(j)-0.05/dtV):ALLpeaktime(j),i))-AHP));
                AHP1pos=(ALLpeaktime(j)-0.05/dtV)+AHP1pos;
                Risepoint=(waveformV(ALLpeaktime(j),i)+waveformV(AHP1pos,i))/2;
                [risepoint,risepos]=min(abs(waveformV(ALLpeaktime(j)-0.05/dtV:ALLpeaktime(j),i)-Risepoint));
                Risepos(j)=risepos+ALLpeaktime(j)-0.05/dtV;
                Decaypoiont=(waveformV(ALLpeaktime(j),i)+waveformV(AHPpos(1,j),i))/2;
                [decaypoint,decaypos]=min(abs(waveformV(ALLpeaktime(j): AHPpos(1,j),i)-Decaypoiont));
                Decaypos(j)=decaypos+ALLpeaktime(j);
                CpeakHW(j)=(Decaypos(j)-Risepos(j))*dtV; 
                
            elseif length(ALLpeaktime)>=2 % If there are two or more peaks
                AHP=[];
                AHPposition=[];
                AHPpos=[];
                AHP1piont=[];
                AHP1pos=[];
                Decaypos=[];
                Risepos=[];
                CpeakHW=[];
                firstHW=[];
                for j=1:length(ALLpeaktime)
                    if j==1
                        [AHP(j),AHPposition]=min(waveformV(ALLpeaktime(j):(ALLpeaktime(j)+0.05/dtV),i));
                        AHPpos(1,j)=AHPposition+ALLpeaktime(j);
                        [AHP1piont,AHP1pos]=min(abs((waveformV((ALLpeaktime(j)-0.03/dtV):ALLpeaktime(j),i))-AHP));
                        AHP1pos=(ALLpeaktime(j)-0.03/dtV)+AHP1pos;
                        Risepoint=(waveformV(ALLpeaktime(j),i)+waveformV(AHP1pos,i))/2;
                        [risepoint,risepos]=min(abs(waveformV(ALLpeaktime(j)-0.03/dtV:ALLpeaktime(j),i)-Risepoint));
                        Risepos(j)=risepos+ALLpeaktime(j)-0.03/dtV;
                        risepos=[];
                        Decaypoiont=(waveformV(ALLpeaktime(j),i)+waveformV(AHPpos(1,j),i))/2;
                        [decaypoint,decaypos]=min(abs(waveformV(ALLpeaktime(j): AHPpos(1,j),i)-Decaypoiont));
                        Decaypos(j)=decaypos+ALLpeaktime(j);
                        decaypos=[];
                        CpeakHW(j)=(Decaypos(j)-Risepos(j))*dtV; 
                        
                    elseif j>=2 & j<=length(ALLpeaktime)-1
                        [AHP(j),AHPposition]=min(waveformV(ALLpeaktime(j):ALLpeaktime(j+1),i));
                        AHPpos(1,j)=AHPposition+ALLpeaktime(j);   
                        Risepoint=(waveformV(ALLpeaktime(j),i)+waveformV(AHPpos(1,j-1),i))/2;
                        [risepoint,risepos]=min(abs(waveformV(AHPpos(1,j-1):ALLpeaktime(j),i)-Risepoint));
                        if ~isempty(risepos)
                            Risepos(j)=risepos+AHPpos(1,j-1);
                            Decaypoiont=(waveformV(ALLpeaktime(j),i)+waveformV(AHPpos(1,j),i))/2;
                            [decaypoint,decaypos]=min(abs(waveformV(ALLpeaktime(j): AHPpos(1,j),i)-Decaypoiont));
                            Decaypos(j)=decaypos+ALLpeaktime(j);
                            CpeakHW(j)=(Decaypos(j)-Risepos(j))*dtV;
                        end
                    else
                        lastime=[];
                        if ALLpeaktime(j)+0.05/dtV>pretime/dtV+0.5/dtV
                          lastime= pretime/dtV+0.5/dtV;
                        else
                          lastime=ALLpeaktime(j)+0.05/dtV; 
                        end
                        [AHP(j),AHPposition]=min(waveformV(ALLpeaktime(j):lastime,i))
                        if ~isempty(AHPposition)
                            AHPpos(1,j)=AHPposition+ALLpeaktime(j);   
                            Risepoint=(waveformV(ALLpeaktime(j),i)+waveformV(AHPpos(1,j-1),i))/2;
                            [risepoint,risepos]=min(abs(waveformV(AHPpos(1,j-1):ALLpeaktime(j),i)-Risepoint));
                            if ~isempty(risepos)
                                Risepos(j)=risepos+AHPpos(1,j-1)
                                Decaypoiont=(waveformV(ALLpeaktime(j),i)+waveformV(AHPpos(1,j),i))/2;
                                [decaypoint,decaypos]=min(abs(waveformV(ALLpeaktime(j): AHPpos(1,j),i)-Decaypoiont));
                                Decaypos(j)=decaypos+ALLpeaktime(j);
                                CpeakHW(j)=(Decaypos(j)-Risepos(j))*dtV;
                            end
                        else
                        end
                    end
                end
            end
          
            slowpeak=find(CpeakHW>=0.005); % Find peaks with wide half-width
            cpeaktime=ALLpeaktime(slowpeak); % Find wide peaks
            AHPpos=AHPpos(slowpeak);
            Risepos=Risepos(slowpeak);
            Decaypos=Decaypos(slowpeak);
            Slowpeak=find(waveformV(cpeaktime,i)-waveformV(AHPpos,i)>20);
            Risepos=Risepos(Slowpeak);
            Decaypos=Decaypos(Slowpeak);
            cpeaktime=cpeaktime(Slowpeak);
            AHPpos=AHPpos(Slowpeak);
            deletedCP=[]
            for s=1:length(speaktime)
               for ss=1:length(cpeaktime)
                 if abs(cpeaktime(ss)-speaktime(s))<10
                     deletedCP=cpeaktime(ss)
                 end
               end
            end
           
            cpeaktime=setdiff(cpeaktime,deletedCP)
            figure,clf
            plot(tspanV,waveformV(:,i)),hold on
            plot(tspanV(cpeaktime),waveformV(cpeaktime,i),'ro')
            % plot(tspanV(AHPpos),waveformV(AHPpos,i),'bo')
            % plot(tspanV(Risepos),waveformV(Risepos,i),'co')
            % plot(tspanV(Decaypos),waveformV(Decaypos,i),'go')
            plot(tspanV(round(speaktime)),waveformV(round(speaktime),i),'ko')
            drawnow
            
            pause
            close

           cpeaktime(cpeaktime*dtV>(Dura(i)+pretime-0.001)| cpeaktime*dtV<(pretime+0.008))=[]% remove trailing stimulus peaks
           CspkFreq(i)=length(cpeaktime)/0.5;
           CspkTime{i}=cpeaktime;
           spkFreq(i)=length(speaktime)/0.5;
           spkTime{i}=speaktime;
           
           % Plot to check identified sodium peaks
           figure,clf
           subplot(3,1,1)
           addone=zeros(1,1)
           SlopwaveformV=[slopwaveformV;addone]
           plot(tspanV,SlopwaveformV),hold on
           title('Slope');
           plot(tspanV(round(sodiumpeaktime)),SlopwaveformV(sodiumpeaktime),'co')
           plot(tspanV(round(calciumpeaktime)),SlopwaveformV(calciumpeaktime),'ro')
           subplot(3,1,2)
           plot(tspanV,waveformV(:,i)),hold on
           title('spike type');
           plot(tspanV(round(sodiumpeaktime)),waveformV(round(sodiumpeaktime),i),'co')
           plot(tspanV(round(speaktime)),waveformV(round(speaktime),i),'ko')
           plot(tspanV(round(calciumpeaktime)),waveformV(round(calciumpeaktime),i),'ro')
           plot(tspanV(round(cpeaktime)),waveformV(round(cpeaktime),i),'go')
           subplot(3,1,3)
           plot(tspanI,waveformI(:,i))
           title('Stimulation'); 
           drawnow
           pause
           close

        end   
        temp=waveformV-repmat(Vrest,size(waveformV,1),1);
        VAHP=trapz(temp(tspanV>0.53&tspanV<0.73,:),1);
       
        % StimAmp
        StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));
        StimAmp(StimAmp<200)=round(StimAmp(StimAmp<200)/10)*10;
        StimAmp(StimAmp>=200)=round(StimAmp(StimAmp>=200)/50)*50;
        APcurrentThrehold = find(spkFreq,1)%20220818
        APcurrentThrehold = StimAmp(1,APcurrentThrehold)%20220818
        % Remove repeated stimuli?
        [StimAmp,x_position]=unique(StimAmp)
        spkFreq= spkFreq(x_position)
        
        CspkFreq=CspkFreq(x_position)
        Adaptation=nan(length(StimAmp),1);
        for i=1:length(StimAmp)
            if length(spkTime{i})>2
                Adaptation(i)=1-(diff(spkTime{i}(1:2)))./diff(spkTime{i}(end-1:end));
            end
        end
        save(['DATA_FIcurve_',filename],'tspanV','tspanI','waveformV','waveformI','TimeStim','Dura',...
          'slopesodium', 'CspkFreq','CspkTime','Vrest','VAHP','spkFreq','spkTime','StimAmp','APcurrentThrehold');
    end
end
