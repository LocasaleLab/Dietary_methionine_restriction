%------------------------------------------------------------------------
%Read and clean up the fold change data
%------------------------------------------------------------------------
Samples={'BC231 plasma','C57BL/6J','Sarcoma','Human',...
    'CRC119','CRC240'};
Samples_TwoLines={sprintf('BC231'),'C57BL/6J',...
    sprintf('Sarcoma'),'Human',...
    sprintf('CRC119'),sprintf('CRC240')};
load FCs.txt
IdxOrder=[4 2 5 6 3];
FCs=FCs(:,IdxOrder);
Samples=Samples(IdxOrder);
Samples_TwoLines=Samples_TwoLines(IdxOrder);
Mets=textread('MetaboliteOverlap.txt','%s','delimiter','\n');
count=0;

%------------------------------------------------------------------------
%Load the lists of methionine-related metabolites with different distance
%cutoffs
%------------------------------------------------------------------------
MS_Met3=textread('MS_Met3.txt','%s','delimiter','\n');
MS_Met4=textread('MS_Met4.txt','%s','delimiter','\n');
MS_Met5=textread('MS_Met5.txt','%s','delimiter','\n');

[a bfc3]=ismember(Mets,MS_Met3);
[a bfc4]=ismember(Mets,MS_Met4);
[a bfc5]=ismember(Mets,MS_Met5);

%------------------------------------------------------------------------
% Draw heatmaps and scatter plots to compare the fold changes between
% models
%------------------------------------------------------------------------
tag=bfc4;
[CovMat,PCovMat]=corr(FCs(tag>0,:),'type','Spearman');
figure;
count=0;
for i=1:4
    for j=i+1:5
        count=count+1;
        subplot(4,5,(j-1)*5+i-5);
        hold on;
        box on;
        scatter(log10(FCs(tag>0,i)),log10(FCs(tag>0,j)),200,'Marker','.','MarkerEdgeColor',[0.85 0.33 0.1]);
        xlim([-1 1]);ylim([-1 1]);
        if j==5
            xlabel(Samples_TwoLines(i));
        end
        if i==1
            ylabel(Samples_TwoLines(j));
        end
        CovVec(count,1)=CovMat(i,j);
        CompareVec(count,1)=strcat(Samples(i),',',Samples(j));
    end
end
colormap('redbluecmap');
figure;
heatmap(CovMat,Samples,Samples,'%0.3f%','MinColorValue',-1,'MaxColorValue',1,...
    'FontSize',14,'TickFontSize',14,'ShowAllTicks',true,'TickAngle',45,'GridLine',':');
colormap('redbluecmap');
title('Spearman Correlation');
figure;
heatmap(PCovMat,Samples,Samples,'%0.1e%','MinColorValue',0,'MaxColorValue',1,...
    'FontSize',14,'TickFontSize',14,'ShowAllTicks',true,'TickAngle',45,'GridLine',':');
colormap('redbluecmap');
title('p-value');