function Io = calc_info(img_seg, groundTruth) %, debug)

% Syntax: I = calc_info(img_seg, truth, debug);
%
%

    
    for j = 1:numel(groundTruth)
    
        truth = groundTruth{1}.Segmentation;
        N = numel(img_seg); % number of pixels in image. (img_seg_Q is log of dominant eigenvector)

        Pa = numel(find(img_seg==1))/N; % probability of a pixel in image being assigned to segment A (by our segmentation).
        Pb = 1-Pa; % probability of pixel being assigned to segment B (total prob = 1) (must be assigned to A or B)
        Is = -N*( Pa*log2(Pa) + Pb*log2(Pb) ); % information contained in our Segmentation of the image.
                                                      % (zero if all pixels in one segment & max if Pa = Pb = 0.5)

        segs = unique(truth); % find all the unique values in the truth matrix
        num_segs = numel(segs); % find number of elements in segs vector (number of segments)

        for k = 1:num_segs % loop through all segments in a single human truth image
            human_seg = find(truth==segs(k));
            size_human_seg(k) = numel(human_seg);
            
            P_ak = numel(find(img_seg(human_seg)==1)) / size_human_seg(k);
            P_bk = 1-P_ak;
            
            Ie(k) = -( P_ak*log2(P_ak) + P_bk*log2(P_bk) ); 
            Ie(isnan(Ie))=0;

        end
         
        
        I(j) = (1/N)*( Is - sum(size_human_seg.*Ie) );  

%         figure, subplot(211),imagesc(img_seg),subplot(212),imagesc(groundTruth{j}.Segmentation)

    end
    
    Io = mean(I);
%     keyboard
    
end

