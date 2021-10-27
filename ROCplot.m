% ROC
realSample=compareAcc(:,1);
possSample=compareAcc(:,2);
ROCfpr=[];
ROCtpr=[];

for threshold=0:0.01:1
    predSample=(possSample>threshold);
    [TPR,FPR]=TPcalcu(realSample,predSample);
    ROCfpr=[ROCfpr FPR];
    ROCtpr=[ROCtpr TPR];
end

figure(3)
plot(ROCfpr,ROCtpr,'b-*');
grid on
hold on
plot([0,1],[0,1],'r-');


function [TPR,FPR]=TPcalcu(realSample,predSample)

trueSheet=realSample+predSample;
TP=sum(trueSheet==2);
TN=sum(trueSheet==0);

falseSheet=realSample-predSample;
FP=sum(falseSheet==-1);
FN=sum(falseSheet==1);

condPositive=TP+FN;
condNegative=FP+TN;

FPR=FP/condNegative;
TPR=TP/condPositive;

end
