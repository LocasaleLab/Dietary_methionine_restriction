%------------------------------------------------------------------------
% Compare fold changes in metabolite levels between different tissue
% types in the sarcoma model
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%Read and clean up the fold change data
%------------------------------------------------------------------------
T_FCs=readtable('FCs_Overlap.txt');
T_Ps=readtable('Ps_Overlap.txt');
Samples=T_FCs.Properties.VariableNames;
Samples=Samples(2:end);
Samples=strrep(strrep(Samples,'_','-'),'PDX','CRC');
Samples_TwoLines=cellfun(@sprintf,strrep(Samples,'_','\n'),'UniformOutput',false);

FCs=table2array(T_FCs(:,2:end));
Ps=table2array(T_Ps(:,2:end));
Mets=table2cell(T_FCs(:,1));

[a idx]=sort(Samples);
Samples=Samples(idx);
FCs=FCs(:,idx);
pPs=-log10(Ps(:,idx));
Samples_TwoLines=Samples_TwoLines(idx);

%-------------------------------------------------------------------------
%Load the lists of methionine-related metabolites with different distance
%cutoffs
%-------------------------------------------------------------------------
MS_Met3=textread('MS_Met3.txt','%s','delimiter','\n');
MS_Met4=textread('MS_Met4.txt','%s','delimiter','\n');
MS_Met5=textread('MS_Met5.txt','%s','delimiter','\n');
count=0;
[a bfc3]=ismember(Mets,MS_Met3);
[a bfc4]=ismember(Mets,MS_Met4);
[a bfc5]=ismember(Mets,MS_Met5);

%-------------------------------------------------------------------------
%Compare fold changes between tissue types
%-------------------------------------------------------------------------
allpos=1:length(Mets);
tag=allpos;%Do the comparison for all metabolites
SubPosList={1:3,4:6};
for npos=1:2
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
            scatter(log10(FCs(tag>0,SubPos(i))),log10(FCs(tag>0,SubPos(j))),'filled','Marker','o');
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

%----------------------------------------------------------------------------------------
%Calculate fration of altered metabolites in methionine related/unrelated metabolites
%----------------------------------------------------------------------------------------
ChangeTags=zeros(size(FCs));
ChangeTags(log(FCs)<0 & pPs>log10(2)+1)=-1;
ChangeTags(log(FCs)>0 & pPs>log10(2)+1)=1;
tag=bfc4;
for i=1:6
    CT=ChangeTags(:,i);
    METN_C(i)=sum(abs(CT(tag>0)));
    METN_NC(i)=sum(tag>0)-METN_C(i);
    Other_C(i)=sum(abs(CT(tag==0)));
    Other_NC(i)=sum(tag==0)-Other_C(i);
    [h,P_FS_right(i),s]=fishertest([METN_C(i) METN_NC(i);Other_C(i) Other_NC(i)],'Tail','right');
    [h,P_FS_left(i),s]=fishertest([METN_C(i) METN_NC(i);Other_C(i) Other_NC(i)],'Tail','left');
	P_FS=min([P_FS_left;P_FS_right]);
    MET_ChangeRatio(i)=METN_C(i)/(METN_C(i)+METN_NC(i)); 
	%Fraction of altered metabolites in methionine-related metabolites
    NoMET_ChangeRatio(i)=Other_C(i)/(Other_C(i)+Other_NC(i));
	%Fraction of altered metabolites in methionine-unrelated metabolites
end