function [numClust, centInd] = decisionGraph(rho, delta, isManualSelect)
%%DECISIONGRAPH Decision graph for choosing the cluster centroids.
%   INPUT:
%       rho: local density [row vector]
%       delta: minimum distance between each point and any other point with higher density [row vector]
%       isManualSelect: 1 denote that all the cluster centroids are selected manually, otherwise 0
%  OUTPUT:
%       numClust: number of clusters
%       centInd:  centroid index vector

    NE = length(rho);
    numClust = 0;
    centInd = zeros(1, NE);
    
    if isManualSelect == 1
        
        fprintf('Manually select a proper rectangle to determine all the cluster centres (use Decision Graph)!\n');
        fprintf('The only points of relatively high *rho* and high  *delta* are the cluster centers!\n');
        figure
        plot(rho, delta, 's', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');
        title('Decision Graph', 'FontSize', 17);
        xlabel('\rho');
        ylabel('\delta');
        
        % original code: selection disappears
        %position = getrect;
        
        % DEW code: highlight selection in red 
        draw_rect = imrect(gca);
        position = wait(draw_rect);
        delete(draw_rect);
        r = rectangle('Position', [position(1), position(2), position(3), position(4)],'Linewidth', 1);
        set(r,'edgecolor','r')
        
        minRho = position(1);
        minDelta = position(2);    
        
        for i = 1 : NE
            if (rho(i) > minRho) && (delta(i) > minDelta)
                numClust = numClust + 1;
                centInd(i) = numClust;
            end
        end
        
    else
        alpha=2;
        P = polyfit(rho,delta,alpha);
        if alpha==1
            figure;
            plot(rho, delta, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
            hold on
            yfit = P(1)*rho+P(2);
            plot(rho,yfit,'g-.');
        elseif alpha==2
            figure
            t2 = 0:1:300;
            y2 = polyval(P,t2);
            plot(rho,delta,'o',t2,y2)
        end
        title('Local-Density Clustering: Decision Graph', 'FontSize', 12);
        xlabel('\rho (Local Density)');
        ylabel('\delta (Minimum Distance)');
        f = polyval(P,rho);
        err=delta-f;
        cutoff=quantile(delta-f,.99);
        for i = 1 : NE
            if err(i)>cutoff
                numClust = numClust + 1;
                centInd(i) = numClust;
            end
        end
    end

end