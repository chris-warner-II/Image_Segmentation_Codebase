function [] = GaussianBoxLoopGen(boxDims,numRuns)

% Create an ENSEMBLE OF IMAGES with an inside and outside box drawn from
% different NORMAL DISTRIBUTIONS.  These images will be used as input to 
% the different segmentation algorithms to compare their performance.


% OTHER IDEAS FOR TESTS::
% INSIDE BOX HAS SAME MEAN BUT MUCH TIGHTER VARIANCE THAN OUTSIDE...
% INSIDE BOX MOVES AROUND INSIDE OUTSIDE BOX TO INVESTIGATE EDGE EFFECTS...
% MORE COMPLICATED EDGES AND CONTOURS THAN JUST BOX... (Curvature effects segmentation?)
    
% Parameters for stochasticity
rng(1234567); % Seed for Random Number Generator

% Check if you are using cluster to adjust home directory
dirPre = onCluster;

for k = 1:numel(boxDims)
    
    ximg = boxDims(k); 
    yimg = boxDims(k); 


    % Make directories to house stacks of input images generated
    if ~exist([dirPre,'images/GaussianBox/',num2str(ximg),'x',num2str(yimg)],'dir')
        mkdir([dirPre,'images/GaussianBox/',num2str(ximg),'x',num2str(yimg)])
    end

    x2 = round(0.4*ximg); % start and end of x&y of inner rectangle
    x3 = round(0.7*ximg); % middle is 0.4 0.7
    y2 = round(0.4*yimg);
    y3 = round(0.7*yimg);

    % Set Means of distributions to be 0 and 1.
    mu=0;
    mu1 = mu;    % Mean of Distribution of inside box
    mu2 = 1-mu;  % Mean of Distribution of outside box

    % Increase variance of distributions to be larger - increase overlap.
    sig = [1e-10, 0.01, 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4];

    % CAN ALSO CHANGE MEAN WHILE HOLDING VARIANCE CONSTANT... FOR DIFFERENT VARIANCES
    % COULD ALSO VARY VARIANCE FOR DIFFERENT MEAN VALUES [0.1 - 0.9 :: 0.2 - 0.8]


    for i = 1:numel(sig)
        sig1 = sig(i);  % Spread of Distribution of inside box
        sig2 = sig(i);  % Spread of Distribution of outside boxlets j


        disp(['Input Image is ',num2str(ximg),'x',num2str(yimg),' box with pixel intensity vals distribution different from surrounding.']);
        iids = struct('name','gaussBox    '); % note: spaces after are important


        for j = 1:numRuns
            imEns(:,:,j) = GaussianBox(ximg,yimg,x2,x3,y2,y3,mu1,sig1,mu2,sig2,0);
        end

        fOut = ['gauss_',num2str(ximg),'x',num2str(yimg),'x',num2str(numRuns),'_Box_',num2str(mu1,'%0.5g'),'_',num2str(sig1,'%0.5g'),'_',num2str(mu2,'%0.5g'),'_',num2str(sig2,'%0.5g')];
        fOut = fOut(fOut~='.');

        gndTruth{1} = GaussianBox(ximg,yimg,x2,x3,y2,y3,0,0,1,0,0); % This creates a binary, black & white, box.
        gndTruth{1}(gndTruth{1}==0) = 2;                            % so regions in ground truth are labeled 1 & 2.

        % Build Distributions of Pixel Intensity Values for inside & outside boxes
        x = linspace(0,1);
        G1 = (1/(sig1.*sqrt(2*pi))).*exp( -(x-mu1).^2 ./ (2.*sig1.^2) );
        G2 = (1/(sig2.*sqrt(2*pi))).*exp( -(x-mu2).^2 ./ (2.*sig2.^2) );
        G1 = G1./sum(G1); G2 = G2./sum(G2); % Normalize (a hack)

        %KL Divergence as measure of difference between distributions
        div = G2./G1;
        div(isnan(div))=1;
        D_KL = sum(log(div).*G2); % = 0 if same; = large if different
        if(isnan(D_KL))
            D_KL = 500; % This is infinity and means no overlap with distributions
        end

        % Plot and save 2 Gaussian Distributions of pixel values to view overlap.
        if(0)  
            hgd = figure; hold on, plot(x,G1,'r','LineWidth',2), plot(x,G2,'b','LineWidth',2)
            title(['Input PIV Distrib: \mu=[',num2str(mu1),',',num2str(mu2),'] \sigma=[',num2str(sig1),',',num2str(sig2),']'],'FontSize',20,'FontWeight','Bold');
            xlabel('Pixel Intensity Value','FontSize',16,'FontWeight','Bold')
            ylabel(['KL Divergence = ',num2str(D_KL)],'FontSize',16,'FontWeight','Bold')
            legend('Inside rect','Outside rect','Location','North')
            saveas(hgd,['./images/GaussianBox/',fOut],'jpg');
            close
        end

        % save each image stack to a mat file to be used as input later.
        save([dirPre,'images/GaussianBox/',num2str(ximg),'x',num2str(yimg),'/',fOut],'imEns','gndTruth','D_KL','mu1','sig1','mu2','sig2','ximg','yimg','x2','x3','y2','y3','numRuns');

    end

end

