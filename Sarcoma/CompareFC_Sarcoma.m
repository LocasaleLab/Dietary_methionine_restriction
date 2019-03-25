%----------------------------------------------------------------------
% Compare fold changes in metabolite levels between different treatment
% conditions in the sarcoma model
%----------------------------------------------------------------------

%----------------------------------------------------------------------
% Read fold changes and sample informations from files
%----------------------------------------------------------------------
Samples=textread('SampleList.txt','%s','delimiter','\n');
load FCs.txt
[a idx]=sort(Samples);
Samples=Samples(idx);
FCs=FCs(:,idx);
Mets=textread('MetaboliteOverlap.txt','%s','delimiter','\n');
count=0;
allpos=1:length(Mets);
tag=allpos;
SubPosList={1:3,4:6,7:9,10:12};

%-----------------------------------------------------------------------
% Draw heatmaps and scatter plots to compare the fold changes between 
% treatment groups
%-----------------------------------------------------------------------
for npos=1:4
    SubPos=SubPosList{npos};
    figure;
    title('Fold Change Comparison');
    count=0;
    for i=1:2
        for j=i+1:3
            count=count+1;
            subplot(1,3,count);
            hold on;
            box on;
            scatter(log10(FCs(tag>0,SubPos(i))),log10(FCs(tag>0,SubPos(j))),'filled',...
                'Marker','o','MarkerEdgeColor',[0.85 0.33 0.1],'MarkerFaceColor',[0.85 0.33 0.1]);
            xlabel(Samples(SubPos(i)));
            ylabel(Samples(SubPos(j)));            
            xlim([-1 1]);
            ylim([-1 1]);
        end
    end
    colormap('redbluecmap');
    CovMat=corr(FCs(tag>0,SubPos),'type','Spearman');
    figure;
    heatmap(CovMat,Samples(SubPos),Samples(SubPos),'%0.3f%','MinColorValue',-1,'MaxColorValue',1,...
        'FontSize',14,'TickFontSize',14,'ShowAllTicks',true,'TickAngle',45,'GridLines',':');
    title('Spearman Correlation');
    colorbar;
    colormap('redbluecmap');
end