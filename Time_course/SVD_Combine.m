%-------------------------------------------------------
% Perform singular value decomposition on the time-course
% data in healthy mice fed control/MR diets
%-------------------------------------------------------

%-------------------------------------------------------
% Read and clean the time course metabolomics data
%-------------------------------------------------------
load Raw_Con.csv
load Raw_MR.csv
logCon=log10(Raw_Con);
logMR=log10(Raw_MR);
logComb=[logCon logMR];

%-------------------------------------------------------
% Calculate singular value decomposition (SVD) of the 
% time-course metabolomics data and plot singular values
%-------------------------------------------------------
XComb=RowNorm(ColNorm(RowNorm(logComb)));
[UComb,SComb,VComb]=svd(XComb);
figure;
plot(diag(SComb),'-o','LineWidth',2);
title('Singular Values of Metabolomics Profile');
xlabel('N');
ylabel('Nth Singular Value');

%--------------------------------------------------------
% Calculate the modes and plot the top three modes with
% highest singular values
%--------------------------------------------------------
Modes_Comb=SComb*VComb';
nkeep=3;
figure;
Time=[0 1 2 4 7 10 14 17 21];
for i=1:nkeep
    subplot(nkeep,1,i);
    XCon=reshape(Modes_Comb(i,1:45),5,9);
    XMR=reshape(Modes_Comb(i,46:end),5,9);
    for j=1:9
        [h PT2Mat(i,j)]=ttest(XCon(:,j),XMR(:,j));
    end
    plot(Time,reshape(Modes_Comb(i,1:45),5,9),'o','LineWidth',2,'MarkerEdgeColor',[0 0 0]);
    hold on;
    h1=plot(Time,mean(reshape(Modes_Comb(i,1:45),5,9)),'LineWidth',2,'Color',[0 0 0]);
    plot(Time,reshape(Modes_Comb(i,46:end),5,9),'o','LineWidth',2,'MarkerEdgeColor',[1 0 0]);
    hold on;
    h2=plot(Time,mean(reshape(Modes_Comb(i,46:end),5,9)),'LineWidth',2,'Color',[1 0 0]);
    title(strcat('Mode ',sprintf('%d',i)));
    legend([h1 h2],{'Control','MR'});
    if i==nkeep
        xlabel('Day');
    end
end

%-------------------------------------------------------
% Functions transforming columns or rows of a matrix to 
% make their mean values zero
%-------------------------------------------------------
function XNorm=ColNorm(X)
	XNorm=X-mean(X);
end

function XNorm=RowNorm(X)
	XNorm=X-mean(X,2);
end