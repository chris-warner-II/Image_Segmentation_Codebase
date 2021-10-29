function [W, Wconn, Wunc, Wdist, Mask] = calc_weights(imagin,sigpix,sigdist,maskFlg,rmax,plt)

% syntax: [W, Wconn, Wunc] = calc_weights(imagin,sigpix,sigdist,threshdist)
%
% This function takes in an input image and some parameters and calculates
% and outputs an Adjancency Matrix for pixels within that image. 
%
% INPUT: (imagin,sigpix,sigdist,threshdist)
%    imagin  - input image
%    sigpix  - sigma value on a gaussian nonlinearity that converts pixel difference abs(Zi - Zj) \in [0,1] into another number \in [0,1]
%    sigdist - sigma value on a gaussian nonlinearity that converts pixel distance abs(i-j) \in [0,1] into another number \in [0,1]
%    maskFlg -
%    distBin -
%    rmax    -
%    plt     - flag to display things and keyboard for error checking.
%
% OUTPUT: [W, Wconn, Wunc]
%   W - Adjacency Matrix using both connectivity constraints and pixel differences
%   Wconn - "Hardware" Connectivity constraints imposed by programmer
%   Wunc - Weights based solely on pixel similarity (unconstrained by distance between pixels)

% preallocate weights matrix (will be mostly zeros because of nearest neighbor connectivity)
W = sparse( numel(imagin) , numel(imagin) );
Wdist = zeros( numel(imagin) , numel(imagin) );
Wpix = zeros( numel(imagin) , numel(imagin) );

%% Calculate Wdist. given user input values for sigdist & threshdist
y_vec = [1:size(imagin,1)]';
y_mat = repmat(y_vec,1,size(imagin,2));
y_unwrapt = reshape(y_mat,1,numel(imagin));
y_un_mat = repmat(y_unwrapt,numel(y_unwrapt),1);
%
x_vec = [1:size(imagin,2)];
x_mat = repmat(x_vec,size(imagin,1),1);
x_unwrapt = reshape(x_mat,1,numel(imagin));
x_un_mat = repmat(x_unwrapt,numel(x_unwrapt),1);
%
x_dist = abs(x_un_mat - x_un_mat');  %
y_dist = abs(y_un_mat - y_un_mat');  %
r = sqrt(x_dist.^2 + y_dist.^2);     % Euclidian or L2 distance between two pixels.
%   
Wdist = exp( -( r - 1) ./ (2*sigdist.^2) ); % Pixel proximity (Gaussian fall-off part)
Wdist = Wdist.*(r<=rmax); % set anything beyond rmax to zero




%% Create a Series of Masking Matrices that Index into weights matrix to get pixels that are separated by a
% given distance in the image (equivalence classes).
switch(maskFlg)
    case(0)
        Mask(1).mask=0;
    case(1)
        rd = unique(r);
        for i=1:numel(rd)
            Mask(i).distance = rd(i);
            Mask(i).mask = (r==rd(i));
        end
        
    case(2)     
        % Consider direction or orientation between pixels
        xd = unique(x_dist);
        yd = unique(y_dist);
        k=0;
        for i=1:numel(xd)
            for j=1:numel(yd)
                k=k+1;
                Mask(k).distance = sqrt(xd(i).^2 + yd(j).^2);
                Mask(k).distX = xd(i);
                Mask(k).distY = yd(j);
                Mask(k).mask = ( (x_dist==xd(i)) & (y_dist==yd(j)) );
            end
        end
end



% Then later in PIV_vs_MaskDist, calculate "Diagonal Means" or Distance Pixel Dependence as: B(1) = mean(mean(W(Mask(1).mask)))


%% Rasterize Image vertically to calculate Wpix matrix (similarity of Pixel Intensities). 
% Turn imagin pixel values into a vector going down column -> imagin(1,1:XXX) [note: same as y dim above].
img_unwrapt = reshape(imagin,1,numel(imagin)); 
img_un_mat1 = repmat(img_unwrapt,numel(img_unwrapt),1);   
img_un_mat2 = repmat(img_unwrapt',1,numel(img_unwrapt));
%                     
Wpix = exp( -( ( (img_un_mat1 - img_un_mat2).^2 ) )./ (2*sigpix^2) ); % Correlation between pixel values

%% Save output Weights Matrix, Hardware Connectivity Constraint Matrix & Pixel Intensity Similarity Matrix.
W = Wdist .* Wpix;
%
W = sparse(W);          % Adjacency Matrix Constrained by Hardware Connectivity
Wconn = sparse(Wdist);  % Hardware Connectivity Component of Adjacency Matrix
Wunc = sparse(Wpix);    % Pixel Similarity Component of Adjacency Matrix


%% Make Plots of Pixel Intensity Dependence on Distance.
if(0)
    
    % Plot Weights (Adjacency) Matrix and the 2 elements that go into its calculation
    figure, hold on
    subplot(131), imagesc(Wpix), title('Unconstrained Pixel Dependence (W_{pix})','FontSize',20,'FontWeight','Bold')
    subplot(132), imagesc(Wdist), title('Hardware Connectivity Constraint (W_{dist})','FontSize',20,'FontWeight','Bold'), axis off
    subplot(133), imagesc(W), title('Constrained Pixel Dependence (W)','FontSize',20,'FontWeight','Bold'), colorbar('East'), axis off


    % make and plot logical matrices of diagonals pertaining to certain distances
    for i = 0:max(size(imagin))-1 % Vertical dimension = direction of image rastering.
        eval(['x',num2str(i),' = (r==',num2str(i),');'])
    end

    figure, imagesc(r),colorbar,title('distance of pixels in image','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x0),colormap('bone'),title('pixels with distance 0','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x1),colormap('bone'),title('pixels with distance 1','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x2),colormap('bone'),title('pixels with distance 2','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x3),colormap('bone'),title('pixels with distance 3','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x0+2*x1+3*x2+4*x3),colormap(bone(4)),colorbar,title('pixels with distance {0,1,2,3}','FontSize',20,'FontWeight','Bold')

    figure, imagesc(x4),colormap('bone'),title('pixels with distance 4','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x5),colormap('bone'),title('pixels with distance 5','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x6),colormap('bone'),title('pixels with distance 6','FontSize',20,'FontWeight','Bold')
    figure, imagesc(x7),colormap('bone'),title('pixels with distance 7','FontSize',20,'FontWeight','Bold')


    % Calculate and plot mean and standard deviation of each diagonal (if STD=0) diagonal is uniform valued
    for d=0:numel(imagin)-1
        dd(d+1,:) = [mean(diag(r,d)),std(diag(r,d)),max(diag(r,d)),min(diag(r,d))];
    %     figure, hold on,
    %     plot(diag(r,d),'g','LineWidth',2); 
    %     line([1, numel(diag(r,d))],[ dd(d+1,1),dd(d+1,1)])
    %     title(['Diagonal #',num2str(d)],'FontSize',20,'FontWeight','Bold'),
    %     xlabel(['Entry #'],'FontSize',18,'FontWeight','Bold'),
    %     ylabel(['Value'],'FontSize',18,'FontWeight','Bold'),
    %     legend('element value','mean')
    end

    figure, hold on,
    title('Diagonal Statistics','FontSize',20,'FontWeight','Bold')
    plot(dd(:,1),'b','LineWidth',2),plot(dd(:,2),'r','LineWidth',2)
    plot(dd(:,3),'g--','LineWidth',2),plot(dd(:,4),'k--','LineWidth',2)
    legend('mean', 'standard deviation','max','min')
    xlabel('Diagonal Number','FontSize',20,'FontWeight','Bold')

    % Find distances in matrix remaining when constant diagonals are removed. Lots of structure left.
    rl = r;
    for i = 0:max(size(imagin))-1 % Vertical dimension = direction of image rastering.
        eval(['rl',num2str(i),'(r==',num2str(i),')=0;'])
    end
    %figure, imagesc(rl),colorbar,title('Remaining Structure in Pixel Distance after Constant Diagonals removed.','FontSize',16,'FontWeight','Bold')


    keyboard

end