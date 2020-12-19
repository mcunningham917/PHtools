%% EstimateRMoScales
% 
% Description
% 
% Use Gaussian Mixture Model (GMM) to find populations of catchments with
% characteristic R_Mo. 
% 
%% Input and output files

supercatchmentFilepath = fullfile(phAnalysisFilePath,groupArea,'Ac0p5km2','PHBs','Cusum02_BenchLength3Steps','Tables')

allPHBs = dir(supercatchmentFilepath);
    j=0;
    
for i = 1:length(allPHBs)
    

    fileName = allPHBs(i).name;
    
    if(~isempty(strfind(fileName, 'Supercatchment')))
        
        j = j+1;
        allPHBList(j) = allPHBs(i);
    end
end

%% Gather list of supercatchment tables to plot from folder

targetList = 1:length(allPHBList)


for  i = 1:length(targetList) 

    targetFolder = targetList(i);
    fileName = allPHBList(targetFolder).name;
    
   
    superWordLocation1 = strfind(fileName, 'Supercatchment');
    superWordLocation2 = strfind(fileName, '_all');
    inBetweenTextSuper = fileName(superWordLocation1+14:superWordLocation2-1);
    superNumFileList(i) = str2double(inBetweenTextSuper); 
    
end

superNumFileList = sort(superNumFileList);
[~,targetList] = find(supercatchmentNum==superNumFileList);

%% Use GMM to estimate relief scales in each supercatchment

for i = 1:length(targetList) 
    clear phbOutletArray modeOutletArray modes outlets R_Mo gmm
    
    superNum = targetList(i);
    
    benchSuperFolderName = ['Supercatchment', num2str(superNum),'PHBs'];
    benchSuperListName = ['Supercatchment', num2str(superNum),'_allOutletModePairs.txt'];
    benchSuperOutName = ['Supercatchment', num2str(superNum)];
    phbOutletArray = dlmread(fullfile(supercatchmentFilepath, benchSuperListName),'\t',2);

    modeOutletArray = phbOutletArray(:,[1,2]);
    modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
    zeroIndices = find(modeOutletArray(:,2)==0);
    modeOutletArray(zeroIndices,:)=[];

    modes = modeOutletArray(:,2);
    outlets = modeOutletArray(:,1);
    R_Mo = modes-outlets;

    
    %Compute GMM for for R_Mo in each supercatchment
    
    %Test for a range of relief scales. 3 are predicted, but more may be
    %present. We choose 10 as a maximum.
    
    numScales=1:10; 
    nK = numel(numScales);
    RegularizationValue = 0.000001;
    options = statset('MaxIter',10000);

        % Fit all models
        for gmmNum = 1:nK;
            gmm{gmmNum} = fitgmdist(R_Mo,numScales(gmmNum),'RegularizationValue',RegularizationValue, 'Options',options);
        end
        
        % Assert  the minimum number of relief scales to be 3, and
        % find whether populations overlap at 1 sigma for each GMM with
        % 3-10 component populations.
        
        clear compPropCond overlapCond...
                nonOverlappingComp minNumComponents_condA minNumComponents_condB...
                possComps
            
        for gmmComponentNum =4:10
   
            numComponents = gmmComponentNum;

            clear props means sigma pdStruct finalComponents componentSort...
                meansSort componentProps sigmasSort 
            

            for j = 1:numComponents
                % Each mean from GMM
                means(j) = gmm{numComponents}.mu(j); 
                % 1 sigma for each mean
                sigma(j) = (sqrt(gmm{numComponents}.Sigma(j)));
                % Component mixing proportion
                mixingProportion(j) = gmm{numComponents}.ComponentProportion(j);   
            end
            
            % Sort results of each GMM by mean (lowest mean to highest)
            [meansSort, componentSortIX] = sort(means,'descend');
            componentProps = mixingProportion(componentSortIX);
            sigmasSort = sigma(componentSortIX);
     
            % Compute +/- 1 sigma for each population
            maximum = meansSort+(sigmasSort); % Vector of each mean + 1s
            minimum = meansSort-(sigmasSort); % Vector of each mean - 1s
            maxMin = [maximum;minimum];
        
            componentSort = sort(componentProps);
            
            %% Conditions for selecting most parsimonious GMM.
        
            % 2) Condition A:  Exclude GMM if any populations overlap at
            % 1s.
            
            % We assume the minimum number of relief scales to be 3. We
            % thus start by checking if the 3 component model has
            % overlapping populations, then iteratively if each higher
            % order GMM does as well.
            
            for m = 1:(gmmComponentNum-1)
                % Vector meansSort is highest-to-lowest means. In the if
                % statement, the first condition checks if each mean-1s is
                % greater than the next highest mean+1s. 
                
                    if((minimum(m) < maximum(m+1))|(maximum(m) < minimum(m+1)))
                        overlapCond(m) = 1;
                    else
                        overlapCond(m) = 0;
                    end
    
                    if(sum(overlapCond)>0)
                        nonOverlappingComp(gmmComponentNum) = (gmmComponentNum-1);
                    else
                        nonOverlappingComp(gmmComponentNum) = 0;
                    end
            end
            
            % B) Condition B: Exclude  GMM if one comp. mixing proportion
            % is <0.05 or if two component populations are both <0.1.
            
            if(componentSort(1)<.05)
                compPropCond(gmmComponentNum) = gmmComponentNum;
            else
                compPropCond(gmmComponentNum) = 0;
            end
            
%             if( componentSort(2)<.1) 
%                 compPropCond(gmmComponentNum) = gmmComponentNum;
%             else
%                 compPropCond(gmmComponentNum) = 0;
%             end
        end
        
        
        nonOverlappingComp = nonOverlappingComp(nonOverlappingComp>0);
        minNumComponents_condA = min(nonOverlappingComp);
        
        compPropCond = compPropCond(compPropCond>0); 
        minNumComponents_condB = (min(compPropCond)-1);
        
        possComps = [minNumComponents_condA, minNumComponents_condB];
        numComponents = min(possComps)
        
        [finalMeans, finalMeansIX] = sort(gmm{numComponents}.mu);
        finalSigmas = sqrt(gmm{numComponents}.Sigma(:));
        finalSigmas = finalSigmas(finalMeansIX);
        [maxFinalMean,maxFinalMeansIX] = max(finalMeans);
        
         %Record all means and sigmas in structure array
    
        GmmOutputStruct(i) = struct('Super', superNum,'NumComp', numComponents,'Means', finalMeans, 'Sigmas', finalSigmas);
    
        %% Plot GMM dist

        clear pdStruct weightedDistStruct
        for n = 1:numComponents
            pd = makedist('Normal','mu', gmm{numComponents}.mu(n), 'sigma', sqrt(gmm{numComponents}.Sigma(n)));
            x = (0:(max(R_Mo)+max(R_Mo)*.5));
            y = pdf(pd, x);
            pdStruct(n) = struct('PDF', y); 
        end
    
        for m = 1:length(pdStruct)
            weightedDistribution = (pdStruct(m).PDF)*(gmm{numComponents}.ComponentProportion(m));
            weightedDistStruct(m) = struct('WeightedPDF', weightedDistribution);
        end
    
        final = vertcat(weightedDistStruct(:).WeightedPDF);
        weightedModelPDF = sum(final(1:m,:));
        
        figure(1)
    	clf;
        hold on;
        set(gca,'FontSize',18, 'fontweight', 'bold')
        plot(weightedModelPDF,x, '-','Color', black, 'Linewidth',4)
        yMax = maxFinalMean + (2.5*finalSigmas(maxFinalMeansIX));
        ylim([0 yMax]);
        xlim([0 .0025]);
 
        
end

%% Plot results
% % Lowest mean assumed to be hillslope scale--light blue.
% % Second lowest mean assumed to be colluvial catchment scale--dark blue.
% % Third mean assumed to be fluvial catchment scale--red.
% % All higher means in yellow.
%    
for m = 1:length(GmmOutputStruct)
    
    numComponents = GmmOutputStruct(m).NumComp;
    
    diseqComponents = numComponents-3;

    figure(2)
    hold on;
    set(gca,'FontSize',18, 'fontweight', 'bold')
    set(gca,'xtick',[0:5:25], 'Linewidth', 1.5)
    set(gca,'ytick',[0:500:2700], 'Linewidth', 1.5)
    ylim([0 2700]);
    xlim([0 25])
    set(gca,'units','pixels')
    x0 = 80;
    y0=60;
    width=700;
    height=400;
    set(gca,'position',[x0,y0,width,height])
    set(gcf,'position',[x0,y0,width+100,height+100])
    upperBound = refline(0, 1500);
    upperBound.Color = 'k';
    lowerBound = refline(0, 1000);
    lowerBound.Color = 'k';

    ylabel('Modal R_{Mo}', 'Fontsize', 25, 'Fontweight', 'bold','Interpreter', 'none')
    xlabel('Supercatchment ID [m]', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'none')


    errorbar(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(1), GmmOutputStruct(m).Sigmas(1),...
        'bo','Color', tiDengSek, 'Linewidth',1)
    plot(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(1),...
        'o','MarkerEdgeColor','black', 'MarkerFaceColor',tiDengSek,'MarkerSize', 10);

    errorbar(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(2), GmmOutputStruct(m).Sigmas(2),...
        'bo','Color', blue, 'Linewidth',1)
    plot(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(2),...
        'o','MarkerEdgeColor','black', 'MarkerFaceColor',blue,'MarkerSize', 10);
    
    errorbar(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(3), GmmOutputStruct(m).Sigmas(3),...
        'bo','Color', redOrange, 'Linewidth',1)
    plot(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(3),...
        'o','MarkerEdgeColor','black', 'MarkerFaceColor',redOrange,'MarkerSize', 10);
        
    if diseqComponents>0;
  
        for i = 1:diseqComponents
            
            n = 3+i;
         
         errorbar(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(n), GmmOutputStruct(m).Sigmas(n),...
            'bo','Color', black, 'Linewidth',1)
        plot(GmmOutputStruct(m).Super,GmmOutputStruct(m).Means(n),...
            'o','MarkerEdgeColor','black', 'MarkerFaceColor',yellow,'MarkerSize', 10);

        end
    end

end

