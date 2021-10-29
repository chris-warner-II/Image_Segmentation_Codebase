function [matchOut,dists] = correspondPixelsB(bmap,gT,maxDist)

%
% OUTPUT:
% match = correspondence between ground truth boundary (gT) and thresholded
%         prob boundary map (bmap) allowing for distances up to maxDist pixels.
%
% INPUTS:
% bmap = unthinned binary thresholded probabalistic boundary map from method.
% gT   = binary boundaries drawn from human ground truthers.
% maxDist = max allowable distance for a correspondence to be found (in units of pixels)
%
% I am writing this function to replace correspondPixels because I do not
% know what it is doing.  It is a black-box mex file that came with the
% benchmark code and its output does not make sense.
%
%
% Note: I save multiple output txt files each time through this function:
%        <image_patch#>_d#_ev2.txt - d denotes slop distance allowed in 
%                                    correspondPixelsB computation.




%% Calculate map of distances from closest gT boundary pixel in gT map
[gTd,gTdI] = bwdist(gT);
% gTd is distance map
% gTdI is index map referring to gT pixel spot.
dists = sort(unique(gTd),'ascend');
xx = find(dists <= maxDist); % dists to loop thru, up to maxDist.
dists = dists(xx);






% plot for visualization & development & error checking
if(0)
    figure
    subplot(131), imagesc(gT), title('gT'), axis square off
    subplot(132), imagesc(bmap), title('bmap'), axis square off
	subplot(133), imagesc(gTd), title('gT dist'), axis square off, colorbar
    disp('ok')
end


% Initialize match with all zeros.
matchOut = zeros(size(bmap,1),size(bmap,2),numel(dists));
match = zeros(size(bmap));



if(0)
    disp(['@ distance of 0 pixel(s), # of edge pixels to go thru: ',num2str(numel(find(gT)))])
end



% First, compute DIRECT PIXEL OVERLAP if any (because it is fastest to compute 
% and the other more complicated computation should be same as direct pixel 
% overlap when there is direct overlap).
dpo = bmap&gT;
match = match + dpo;
matchOut(:,:,1) = match;
gT = gT - dpo;
bmap = bmap - dpo; % this is done with ahit in loops down below.


% plot for visualization & development & error checking 
if(0)
    figure
    subplot(131), imagesc(match),title('match found'), axis square off
	subplot(132), imagesc(gT),title('gT left (D=0)'), axis square off
	subplot(133), imagesc(bmap),title('bmap left'), axis square off
    disp('ok')
end








% Loop thru all distances ( >0 but <=maxDist) .
for i = 2:numel(dists)
    
    % make a distance map for all pixels in image of distance from gT pixels remaining to be accounted for.
    [gTd,gTdI] = bwdist(gT);
    
    
    % for all remaining pixels in gT, move them out by one step to see if
    % now there is some overlap with remaining bmap.
    yy = ( (gTd <=dists(i)) & (gTd > dists(i-1)) ); % Check for Overlap between a Distance Pair
    
    
    zz = single(gTdI).*yy; % this is the remaining boundary pixels of gT moved out by a step (bracketed between 2 distances)
                           % and each one is labeled with the location of the pixel in the original gT if originated from
    
    % VIP NOTE: Have to be careful not to "anihilate" two pixels as we look
    % out at some non-zero distance.  Because each single pixel became two
    % pixels! (use this to avoid that somehow)
    if(0)
        figure, imagesc(zz), colorbar
        title({['GT boundaries moved ',num2str(dists(i)),' pixels out'],...
               ['Entry refers to pixel in original gT that this pixel drifted from']})
        disp('ok')   
    end
    
    indPix = sort(unique(zz));
    indPix = indPix(indPix>0); % sorted unique index of entry in gT matrix that corresponds
    
    if(0)
        disp(['@ distance of ',num2str(dists(i)),'pixel(s), # of edge pixels to go thru: ',num2str(numel(indPix))])
    end
    
    % Loops thru each boundary pixel in the remaining gT and checks
    % for overlap between bmap pixel and gTs slop shifted pixels
    for j = 1:numel(indPix)
    
        pc = find(zz==indPix(j)); % pc = pixel correspondence (all pixels in this gTd matrix a distance d from one pixel in gT matrix)
        
        if ~gT(indPix(j)) % note this should be 1. (gT boundary pixel should still exist)
            continue
        end
            
            
        ahit = find( ~match(pc) & bmap(pc) ); % need match to be zero (that bmap pixel is not already matching with another gT boundary)
                                              % and need bmap to be one (a boundary exists at this pixel correspondence as gT pixels expand out distance d)
        
        % Note: I am only doing this for 1 pixel in ahit (this is not optimized for pathological cases tho)
        if ~isempty(ahit) % an overlap occurred between shifted gTd pixel and bmap boundary
            
            if (numel(ahit)>1)
                
                if(0)
                    disp('Multiple Hits within dist range.  Which to choose? This is the optimization.')
                end
                
                % OK, so here I want to plot pixel option along with remaining gT & bmap pixels.
                if(0)
                    figure
                    subplot(131), imagesc(gT) , title('gT remaining')
                    hold on, 
                    scatter( floor(pc(ahit)./size(gT,1))+1, rem(pc(ahit),size(gT,1)), 'wx')   % moved gT edge & possible pixel correspondence
                    scatter( floor(indPix(j)./size(gT,1))+1, rem(indPix(j),size(gT,1)), 'wo') % original gT edge where possible pc's originated from.
                    axis square
                    set(gca,'Xtick',[],'Ytick',[])
                    %
                    subplot(133), imagesc(bmap), title('bmap remaining')
                    hold on, 
                    scatter( floor(pc(ahit)./size(gT,1))+1, rem(pc(ahit),size(gT,1)), 'wx')   % moved gT edge & possible pixel correspondence
                    scatter( floor(indPix(j)./size(gT,1))+1, rem(indPix(j),size(gT,1)), 'wo') % original gT edge where possible pc's originated from.
                    axis square
                    set(gca,'Xtick',[],'Ytick',[])
                    %
                    subplot(132), imagesc(match), title('match currently')
                    hold on, 
                    scatter( floor(pc(ahit)./size(gT,1))+1, rem(pc(ahit),size(gT,1)), 'wx')   % moved gT edge & possible pixel correspondence
                    scatter( floor(indPix(j)./size(gT,1))+1, rem(indPix(j),size(gT,1)), 'wo') % original gT edge where possible pc's originated from.
                    axis square
                    set(gca,'Xtick',[],'Ytick',[])                    
                end
                
                % Idea: As a "tie-breaker" or to optimize the way we compute distance-dependent pixel correspondence:
                % compute gTd &gTdI with the current pixel (white circle) anihilated and see what other remaining gT pixels 
                % the moved gT pixels or the current pixel (white x's) are closest to and use the moved pixel option (white x) 
                % that is furthest from all remaining other gT boundary pixels.
                gTminus = gT;
                gTminus(indPix(j))=0;
                [gTd_minus,gTdI_minus] = bwdist(gTminus);
                best2take = find( gTd_minus(pc(ahit)) == max(gTd_minus(pc(ahit))) ); % index into which of the white x's is the max distance 
                % away from any other remaining groundTruth edges and thus least likely to take away a possible correspondence from them.
                % note: if distances are same, best2take will have multiple entries and now we pick first one arbitrarily.
                
                if(0)
                    disp(['Multiple matching distances to bmap pixels'])
                    gTd_minus(pc(ahit))
                end
                
                
                if(numel(best2take)>1) % if several choices tie as pixel match.
                
                    % 2nd Idea: As a 2nd layer to the tie-breaker, we can look at distances from current match pixels and choose the
                    % white x that is closest to current match pixels  (to enforce some loose line continuation) - nothing here enforces 
                    % any straightness tho. But tends to keep boundary on one side of original boundary instead of hopping back and forth randomly.
                    [match_d,match_dI] = bwdist(match);
                    best2takeB = find( match_d(pc(ahit(best2take))) == min(match_d(pc(ahit(best2take)))) );

                    %match_d(pc(ahit))
                    
                    if(0)
                        disp(['Distance to closest current match pixel'])
                        match_d(pc(ahit(best2takeB)))
                    end
                    
                    % Need to incorporate best2takeB into best2take
                    % Just take 1st one. No further checks.
                    best2take = best2take(best2takeB(1));
                    
                
                else
                
                	% No further work needed if only single best pixel found.
                    
                end
                
                % No further checks: (1). furthest from remaining bmap pixels. 
                %                    (2). Closest to currently existing match pixels.
                match(pc(ahit(best2take(1)))) = 1; % set the first entry of ahit to be a match
                gT(indPix(j))=0;        % clear the original (unmoved) pixel from gT
                bmap(pc(ahit(best2take(1)))) = 0;  % clear the moved pixel from bmap
                
                
                
                
            else
                
                match(pc(ahit(1))) = 1; % set the first entry of ahit to be a match
                gT(indPix(j))=0;        % clear the original (unmoved) pixel from gT
                bmap(pc(ahit(1))) = 0;  % clear the moved pixel from bmap
            
            end
        end
        
        if(0)
            disp(['ok - j = ',num2str(j),' / ',num2str(numel(indPix))])
        end
        
    end % loop over indPix (the number of pixels in gT boundary still existing.
    
    
    
    % plot for visualization & development & error checking
    if(0)
        figure
        subplot(131), imagesc(match),title('match found'), axis square off
        subplot(132), imagesc(gT),title(['gT left (D=',num2str(dists(i),2),')']), axis square off
        subplot(133), imagesc(bmap),title('bmap left'), axis square off
        disp('ok')
    end
    
    matchOut(:,:,i) = match;

end % loop over (allowed or slop) distances in correspondPixels computation.


if(0)
    disp('ok - Go thru once for each GT.')
end


