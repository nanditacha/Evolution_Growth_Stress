%Initial Distribution
mu=500; %trehalose mean
sigma=1; %trehalose 
nsamples=1000; %size of population
tr0 = mu + randn(nsamples,1) * sigma; %initial distribution of trehalose

Nparam=2;
Nreplicates=50;
Ncycles=500;

%trmut=6; %phenotypic effect of mutations
X = linspace(-7,7,1000);
mum = 0;
sigmam = 3;
skewm = 0;
kurtosism = 700;

mutrate=0.8;

% tr1large=NaN(Nreplicates,Ncycles,nsamples); %after growth
% trqlarge=NaN(Nreplicates,Ncycles,nsamples); %after quiescence
% trsurvlarge=NaN(Nreplicates,Ncycles,nsamples); %after survival
% trlarge=NaN(Nreplicates,Ncycles,nsamples); %after dilution

meang=NaN(Nparam,Nreplicates,Ncycles,2);
meanq=NaN(Nparam,Nreplicates,Ncycles,2);
means=NaN(Nparam,Nreplicates,Ncycles,2);
meand=NaN(Nparam, Nreplicates,Ncycles,2);

popsizelarge=zeros(Nparam,Nreplicates,Ncycles,1);
quilarge=zeros(Nparam,Nreplicates,Ncycles,1);
survivallarge=zeros(Nparam,Nreplicates,Ncycles,1);

lagtimemean=zeros(Nreplicates,Ncycles,2);

growthtime=zeros(Nparam, Nreplicates, Ncycles,1); %to record how long they can grow

%Lag Determination %lagtime=tl0*(1-tr0/trl)
tl0=30; %lag time for no trehalose
trl=1000;%trehalose at which no lag

%growth rate
%growthr=g0*(1-tr/trg)^2
trg=1200; %trehalose for no growth
g0=40;
r0=10^5;

%Quiescence
trq0=400; %trehalose for quiescence probability=50%
gamma_q=100;

%Survival parameters 
trs=850; %trehalose for survival probability=50%
gamma_s=50;

%Dilution
d=30;

for param=1:Nparam %different values of chosen parameter
    g0=20*param;
    disp('g0 is')
    disp(g0)
    for rep=1:Nreplicates %each replicate for given parameter value
        disp('replicate');
        disp(rep);
    
        tr0 = mu + randn(nsamples,1) * sigma; %initial distribution of trehalose
    
        for N=1:Ncycles %different cycles of FT        
            disp("Cycle No")
            disp(N);
            %disp("Growth");
        
            %Lag time determined by trehalose
            lagtime=tl0*(1-tr0/trl); 
            lagtimemean(N,1)=mean(lagtime);
            lagtimemean(N,2)=std(lagtime);
              
            %growth dynamics: limited resource
         
            r=r0;%10^5; %resource
            trtemp=tr0;%list of lag time cells
            trtemp2=[];%list of dividing cells
            sizelag=size(trtemp);
            tr1=[]; %offspring in each generation which are dividing
        
        
            for t=1:50 %growth cycles (50 is only the maximum value of t)
                %disp(t)
                sizelag=size(trtemp,2);
        
                %check which ones start dividing
                indexlist=[];
                if sizelag>0
                    for i=1:sizelag %for all the individuals in the lag time population
                        if lagtime(i)<t %if lagtime is over
                            trtemp2=[trtemp2, trtemp(i)]; %add to list of dividing cells
                            indexlist=[indexlist,i]; %keep track of indices to be deleted
                        end 
                    end
                end
        
                sizelag=size(trtemp,2);
                sizediv=size(trtemp2,2);
           
                for i=1:sizediv %for everyone in growth
                    growthr=ceil(g0*(1-trtemp2(i)/trg)*(1-trtemp2(i)/trg));
                    if growthr==1
                        tr1=[tr1,trtemp2(i)];
                        r=r-1;
                    else
                        for j=1:growthr
                            r=r-growthr;
                            mutrand=rand;
                            if mutrand<mutrate
                                tr1=[tr1,trtemp2(i)+pearsrnd(mum,sigmam,skewm,kurtosism,1,1)];%tr1=[tr1,trtemp2(i)+trmut*randn]; %add two offspring to tr1
                                
                            else
                                tr1=[tr1,trtemp2(i)];%tr1=[tr1,trtemp2(i)+trmut*randn]; %add two offspring to tr1
                            
                            end
                            
                        end

                    end
                    
                end
                
                trtemp(indexlist)=[]; %delete from list of lag time cells
                lagtime(indexlist)=[]; %delete from list of lag times
        
                %for everyone in lag
                r=r-sizelag*1;
                %disp(r)                   
                if r<1000
                    %disp(t);
                    growthtime(param, rep, N,1)=t;
                    break;
                    
                end
                trtemp2=tr1; %all dividing cells in next generation  
        
                %disp(size(trtemp2))
            end
        
            meang(param,rep,N,1)=mean(tr1);
            meang(param,rep,N,2)=std(tr1);
        
            popsize=size(tr1,2);
            popsizelarge(param,rep,N)=popsize;
                
            %Selecting Quiescent Cells: According to Trehalose levels
            trq=[];
            for i=1:popsize
                 quprob=1/(1+exp(-(tr1(i)-trq0)/gamma_q));%1-exp(-tr1(i)/trq0);
                x=rand;
                if quprob>x
                    trq=[trq,tr1(i)];
                end
            end
      
            meanq(param,rep,N,1)=mean(trq);
            meanq(param,rep,N,2)=std(trq);
        
            quilarge(param,rep,N,:)=size(trq,2);
                    
            trsurv=[];%matrix of survivors
            
            for i=1:size(trq,2)
                survpr=1/(1+exp(-(trq(i)-trs)/gamma_s));
                x=rand;
                if survpr>x %add individual to survived list only if survival prob is larger than random number
                    trsurv=[trsurv,trq(i)];
                end      
            end
        
            %disp("Survival");
            %disp(size(trsurv,2));
        
            means(param,rep,N,1)=mean(trsurv);
            means(param,rep,N,2)=std(trsurv);
        
            % disp(mean(trsurv));
            survivallarge(param,rep,N,:)=size(trsurv,2);
        
            tr2=randsample(trsurv,ceil(size(trsurv,2)/d)); %re-inoculation into fresh medium
        
            meand(param,rep,N,1)=mean(tr2);
            meand(param,rep,N,2)=std(tr2);
        
            % disp("Reinoculated population");
            % disp(size(tr2,2));
        
            tr0=tr2;
        
            if mod(N,10)==0
                save("trsim_singleplot.mat");
            end
        end
    end
end

    
save("trsim_singleplot.mat")