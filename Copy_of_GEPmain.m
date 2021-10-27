%GEP programe main v0.9.9
%TEST
%this program is aimed to get regression function via GEP algorithm
%attention:this program apply basic GEP and GEP-RNC const method
%main.m include 3 bodies, datainput model, processing model, conduct model
%decode, select, mutate, IS, cross, 5 parts covered

%1.body established, some bug corrected
%2.select,mutate,IS,RIS,3 cross have added
%3.plot function added

%created by shiro_ryou in 2020/4/14/10/28
%last edited by shiro_ryou in 2020/9/11/21/12

clc
clear
close all
tic

%set debug config
%dbstop if error

%initialize the program
popSize=100;%even#,100 is recomand
iterationMax=170;%500
indexFinish=0.0001;
countFinish=0;

%input gene setting info
geneHead=9;
chromNum=3;%each chrom has 3 genes using '+' to connect
Func=['/' '+' '-' '*' 'q' 'p'];
Fnary=[2 2 2 2 1 1];
Tail=['AB?'];%if variable changes,plz change GEPfitness.m, max is P(16)
%Const=[3.14 2.718 1.414 1.732 2.236 0.0001 0.001 0.01 0.05];%num of C is equal to genetail and <9
%Const=5*(rand(1,9)-0.5);
Const=[1.0964  0.4390   -2.1498    0.5188   -1.3976   -0.4095   -0.8592    1.5455   -0.9755];

geneTail=(max(Fnary)-1)*geneHead+1;%genecost = genetail

pmutate=0.144;pis=0.1;pris=0.1;
pcrosss=0.9;pcrossd=0.4;pcrossg=0.2;
isLength=3;

%input data
load('trangle_100.mat')
try
    load('eliteLib.mat')%这里使用常态化精英基因库，无需再初始化
catch
    eliteLib=[];
end
sourceData=originData';
testData=sourceData(91:110,:);
sourceData=sourceData(1:90,:);

%MEDV is y & others are xi
%plz change setting info in the GEPfitness.m when variable number changes

%body
%initialize parameter
maxfitness=0;
bestfitness=0;
bestchrom='shiro_ryou';
bestindividual='shiro_ryou is the best';
compareAcc=1919810;

%initialize plot function
maxfitPlot=zeros(1,iterationMax);
maxvarPlot=zeros(1,iterationMax);
popfitPlot=zeros(1,iterationMax);

%creat newpop
newpop=GEPnewpop(Func,Tail,Const,geneHead,geneTail,popSize,chromNum);

%do the loop
i=1;
while i<=iterationMax
    
    %caculate fitness
    [fitnessList,varList,maxfitness,maxMathexp,maxchrom,compareAcc]=...
        GEPfitness(newpop,geneHead,geneTail,chromNum,Func,Tail,Fnary,Const,sourceData);
    
    %continue?
    if bestfitness==1000 %未知问题可以选择大于900或更低的数值
        popfitPlot(i)=maxfitness;
        maxfitPlot(i)=bestfitness;
        maxvarPlot(i)=bestvariance;
        bestcompare=compareAcc;
        disp('GEP find the bestchorm(1000), so end loop ahead of schedule')
        break;
    end
    
    if abs(bestfitness-maxfitness)<indexFinish
        countFinish=countFinish+1;
        if countFinish==30
            %plot info
            popfitPlot(i)=maxfitness;
            maxfitPlot(i)=bestfitness;
            maxvarPlot(i)=bestvariance;
            bestcompare=compareAcc;
            disp('GEP reach converge, so end loop ahead of schedule')
            break;
        end
    else
        countFinish=0;
    end
    
    popfitPlot(i)=maxfitness;%plot info
    
    %update the bestindividual
    if maxfitness>bestfitness
        bestfitness=maxfitness;
        bestvariance=1000/bestfitness-1;
        bestindividual=maxMathexp;
        bestchrom=maxchrom;
        bestcompare=compareAcc;
    end
    
    %update iretation plot
    maxfitPlot(i)=bestfitness;
    maxvarPlot(i)=bestvariance;
    
    %dynamic plot test
    cla;
    hold on;
    plot(1:i,popfitPlot(1:i),'b');
    plot(1:i,bestfitness*ones(1,i),'r');
    if i==1
        title('Plot 1: maxfitness of each generation')
        xlabel('genaration')
        ylabel('fitness')
        grid on
    end
    %pause(0.01)
    drawnow
    
    %select individuals that join GEPopration
    newpop=GEPselect(fitnessList,newpop);
    
    %mutate
    newpop=GEPmutate(newpop,geneHead,geneTail,chromNum,Func,Tail,Const,pmutate);
    
    %IS & RIS
    newpop=GEPdis(newpop,geneHead,geneTail,chromNum,isLength,Func,pis,pris);
    
    %cross (single double & gene)
    newpop=GEPtcross(newpop,chromNum,pcrosss,pcrossd,pcrossg);
    
    %HGA strategy
    if ~mod(i,20)
        eliteLib=[eliteLib;bestchrom];
    end
    
    if ~mod(i,40)
        newpop=GEPelite(newpop,eliteLib);
    end
    
    i=i+1;%loop i
end

%result and test plot info
Price=inline(vectorize(bestindividual));
[dataNum,varNum]=size(sourceData);

for i=1:varNum-1
    eval([char(64+i),'=testData','(:,',num2str(i),')',';']);
end
yTest=Price(A,B);
figure(2)
scatter(1:length(yTest),yTest,15,'r','filled')
hold on
scatter(1:length(yTest),testData(:,3),18,'b')
grid on

%保存基因库
save('eliteLib.mat','eliteLib')

toc