function [TH_best, TH_loop, modu] = threshModuLineSearch(xxbig,W,img,N)

% syntax [TH_best, THs] = threshModuLineSearch(xxbig,W,N);
%
% This function will go through a N thresholds (including 0 and the mean 
% pixel value) and try segmenting the eigenvector at those values.  It will
% then measure the scalar quality measure, the modularity.  It will output
% the threshold that yielded the largest modularity value, that should lead
% to the best (most meaningful) segmentation.

THs = linspace(min(xxbig(:)),max(xxbig(:)),N);
THs = THs(2:end-1);
TH_loop = [0, mean(xxbig(:)), THs];

for i = 1:numel(TH_loop) % % loop through different values of threshold for Eigenvector

    % find segmentation using threshold on eigenvector
    TH=TH_loop(i);
    seg = xxbig>TH;

    % Calc modularity scalar quality measure
    % Diagonal Wmat Means : calculate modularity scalar (quality of segmentation measure) 
    seg_unwrapt = reshape(seg,1,numel(img)); % unwrapping or rasterizing down the vertical dimension
    seg_unwrapt = 2.*(seg_unwrapt-1/2); % convert SDO to from {0,1} to {-1,+1} so that both segments contribute to modularity value   
    S = sparse(double(seg_unwrapt)'*double(seg_unwrapt));
    S = (S+1)./2; 
    d = sum(W);
    NMng = (d'*d)./sum(d);
    Qng = W - NMng;
    xx=S.*Qng;
    modu(i) = sum(xx(:))./sum(W(:)); 

    if(0)
        % Plot Segmentation for all thresholds for troubleshooting
        segA = img.*seg;
        segB = img.*~seg; % pixel overlays using segments as mask
        %
        h=figure; colormap('bone') % 0
        subplot(223), imagesc(segA),title(['Segment A @ ',num2str(TH)],'FontSize',20,'FontWeight','Bold')
        subplot(224), imagesc( segB ),title(['Segment B modu = ',num2str(modu(i))],'FontSize',20,'FontWeight','Bold')
        subplot(221), imagesc( img ),title(['Image'],'FontSize',20,'FontWeight','Bold')
        subplot(222), imagesc( real(log(xxbig)) ),colorbar,title(['Log E-vec'],'FontSize',20,'FontWeight','Bold')
        set(h,'Position',[1 1 1200 500]);
    end



end


TH_best = TH_loop(modu==max(modu)); % best threshold

TH_best = TH_best(1); % in case >1 threshold got same modularity score