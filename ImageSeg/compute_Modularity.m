function [Q,iter] = compute_Modularity(W,Wconn,Wunc,Mask,ximg,yimg,normlze,maxEnt,maxIter,topo,maskFlg)

% syntax: [Q] = compute_Modularity(W,W_neut,ximg,yimg,lambda,normlze,neg);
%
%   Q is the output Modularity Matrix
%   W is the weights matrix
%   W_conn is a weights matrix generated with uniform input image.  Tell
%     about the connectivity of the network or maximum similarity possible.
%   ximg = horizontal (x) dimension of image - used in null model calculation
%   yimg = vertical (y) dimension of image - used in null model calculation
%   normlze = flag {0,1} to use Modularity normalized by degree of vertices
%   maxEnt = flag {0,1} to use Hillar's MaxEnt algorithm to compute null
%                       model instead of di*dj*bij approximation.


%% Compute Null Model
if(maxEnt)
%     disp('Computing Null Model MaxEnt')
    if ~exist('maxIter','var')
        maxIter = 10000;
    end
    MEth = 1e-6;
    MEdyn =0;
    [NM,iter] = compute_NM_MaxEnt(W,Wconn,ximg,yimg,maxIter,MEth,MEdyn,topo);
else
%     disp('Computing Null Model Heuristic')
    [NM] = compute_null_model(W,Wconn,Mask,topo,maskFlg);
    iter=1;
end

% Compute Modularity from Weights and Null Model
Q = W - NM;








% Maybe look at Normalized Modularity
if(normlze)
    % DO IT HERE IF IT EXISTS! (NormModularity = Q/D or modularity/Degree?)
end


% Some Analysis (Looking at Weights, Null Model & Modularity Matrices vs Image)
if(0)
    figure, imagesc(W), colorbar, title('A')
    figure, imagesc(NM), colorbar, title('N')
    figure, imagesc(Q), colorbar, title('Q')

    Qf = full(Q);
    Wf = full(W);
    Nf = full(NM);

    Wfn = Wf.*(Qf<0)
    Nfn = Nf.*(Qf<0)

    figure,
    subplot(221), imagesc(Qf.*(Qf>0)), colorbar, title('Q>0')
    subplot(222), imagesc(Qf.*(Qf<0)), colorbar, title('Q<0')
    subplot(223), imagesc(Nfn, [min(min(Nfn(:)),min(Wfn(:))),max(max(Nfn(:)),max(Wfn(:)))]), colorbar, title('N(Q<0)')
    subplot(224), imagesc(Wfn, [min(min(Nfn(:)),min(Wfn(:))),max(max(Nfn(:)),max(Wfn(:)))]), colorbar, title('W(Q<0)')

    keyboard
end


%% DONT DELETE - IMPORTANT ANALYSIS BELOW PERHAPS TO FORMALIZE.   
%% Code to compare Topographic ME vs Nontopographic ME vs our Topographic B_{ij} and the original N&G Null Models
if(0)
    
    % Note: I am no longer generating Wunc (the weights unconstrained by distance) so some of this analysis may not work anymore.
    
    % (A). Calculate different topographic and non- Null Models using Topographic Weights Matrix
    NMtme = compute_NM_MaxEnt(W,Wconn,ximg,yimg,maxIter,MEth,MEdyn,1);               % Topographic MaxEnt Null Model (calc with topographic weights)
    NMtme = NMtme.*Wconn;                                            
    NMtng = compute_null_model(W,Wconn,Mask,1);                                      % Topographic New&Giv Null Model (calc with topographic weights)
    NMme = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,maxIter,MEth,MEdyn,0);        % Nontopographic MaxEnt Null Model (calc with topographic weights)
    NMng = compute_null_model(W,ones(size(W)),Mask,0);                               % Nontopographic New&Giv Null Model (calc with topographic weights)



    % (B). Calculate different topographic and non- Null Models using Nontopographic (Unconstrained) Weights Matrix
    NMtme_unc = compute_NM_MaxEnt(Wunc,Wconn,ximg,yimg,maxIter,MEth,MEdyn,1);         % Topographic MaxEnt Null Model (calc with topographic weights)
    NMtme_unc = NMtme.*Wconn;                                            
    NMtng_unc = compute_null_model(Wunc,Wconn,Mask,1);                                % Topographic New&Giv Null Model (calc with topographic weights)
    NMme_unc = compute_NM_MaxEnt(Wunc,ones(size(W)),ximg,yimg,maxIter,MEth,MEdyn,0);  % Nontopographic MaxEnt Null Model (calc with topographic weights)
    NMng_unc = compute_null_model(Wunc,ones(size(W)),Mask,0);                         % Nontopographic New&Giv Null Model (calc with topographic weights)



    %%
    % (A). Plot null models and Topographical Weights along with differnece between null models.
    figure, CHIGH = max([max(NMtme(:)), max(NMtng(:)), max(NMme(:)), max(NMng(:)), max(W(:))]);
    subplot(241), imagesc(NMtme,[0,max(NMtme(:))]); title('Topographic MaxEnt NM','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(242), imagesc(NMme,[0,max(NMme(:))]); title('Nontopographic MaxEnt NM','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(243), imagesc(NMtng,[0,max(NMtng(:))]); title('Topographic N&G NM','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(244), imagesc(NMng,[0,max(NMng(:))]); title('Nontopographic N&G NM','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(245), imagesc(abs(W-NMtme),[0,CHIGH]); title('Adjacency Matrix - Topo MaxEnt','Fontsize',16,'Fontweight','Bold'), colorbar('East'),
    subplot(246), imagesc(abs(W-NMme),[0,CHIGH]); title('Adjacency Matrix - Nontopo MaxEnt','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(247), imagesc(abs(W-NMtng),[0,CHIGH]); title('Adjacency Matrix - Topo N&G ','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(248), imagesc(abs(W-NMng),[0,CHIGH]); title('Adjacency Matrix - Topo MaxEnt','Fontsize',16,'Fontweight','Bold'),colorbar('East'), axis off


    % (B). Plot null models and Unconstrained Weights along with differnece between null models.
    figure, CHIGH = max([max(NMtme_unc(:)), max(NMtng_unc(:)), max(NMme_unc(:)), max(NMng_unc(:)), max(Wunc(:))]);
    subplot(241), imagesc(NMtme_unc,[0,max(NMtme_unc(:))]); title('Topographic MaxEnt NM','Fontsize',16,'Fontweight','Bold'),colorbar('East'), axis off
    subplot(242), imagesc(NMme_unc,[0,max(NMme_unc(:))]); title('Nontopographic MaxEnt NM','Fontsize',16,'Fontweight','Bold'),colorbar('East'), axis off
    subplot(243), imagesc(NMtng_unc,[0,max(NMtng_unc(:))]); title('Topographic N&G NM','Fontsize',16,'Fontweight','Bold'),colorbar('East'), axis off
    subplot(244), imagesc(NMng_unc,[0,max(NMng_unc(:))]); title('Nontopographic N&G NM','Fontsize',16,'Fontweight','Bold'),colorbar('East'), axis off
    %,[0,CHIGH]);
    subplot(245), imagesc(abs(Wunc-NMtme_unc)); title('Adjacency Matrix - Topo MaxEnt','Fontsize',16,'Fontweight','Bold'), colorbar('East'),
    subplot(246), imagesc(abs(Wunc-NMme_unc)); title('Adjacency Matrix - Nontopo MaxEnt','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(247), imagesc(abs(Wunc-NMtng_unc)); title('Adjacency Matrix - Topo N&G ','Fontsize',16,'Fontweight','Bold'), colorbar('East'), axis off
    subplot(248), imagesc(abs(Wunc-NMng_unc)); title('Adjacency Matrix - Topo MaxEnt','Fontsize',16,'Fontweight','Bold'),colorbar('East'), axis off


    %%
    % (A). Plot node incidences (row sums) for MaxEnt vs Heuristic with weights rowsums as correct answer
    figure, hold on, plot(sum(NMtme),'b'),plot(sum(NMme),'r'),plot(sum(NMtng),'g'),plot(sum(NMng),'m'),plot(sum(W),'k')
    title('Nearest Neighbor Adjacency Row Sums (Node Incidences)','Fontsize',16,'Fontweight','Bold')
    legend('Topographic MaxEnt NM','NonTopographic MaxEnt NM','Topographic N&G NM','Nontopographic N&G NM','Weights Mat')
    text(25,5,'Nontopo N&G (mag) on top of Weights (blk)')
    xlabel('Row # of Constrained Weights Matrix','Fontsize',16,'Fontweight','Bold')
    ylabel('Sum of Row','Fontsize',16,'Fontweight','Bold')
    %
    figure, hold on, plot(sum(NMtme),'b'),plot(sum(NMtng),'g'),plot(sum(NMng),'m'),plot(sum(W),'k')
    title('Unconstrained AdjacencyRow Sums (Node Incidences) Zoom [exclude Nontopographic MaxEnt NM]','Fontsize',16,'Fontweight','Bold')
    legend('Topographic MaxEnt NM','Topographic N&G NM','Nontopographic N&G NM','Weights Mat')
    xlabel('Pixel Distance* (* distance in Weights Matrix is distance/orientation combination in image)','Fontsize',16,'Fontweight','Bold')
    ylabel('Expected Weight','Fontsize',16,'Fontweight','Bold')


    % (B). Plot node incidences (row sums) for MaxEnt vs Heuristic with weights rowsums as correct answer
    figure, hold on, plot(sum(NMtme_unc),'b'),plot(sum(NMme_unc),'r'),plot(sum(NMtng_unc),'g'),plot(sum(NMng_unc),'m'),plot(sum(Wunc),'k')
    title('Unconstrained Adjacency Row Sums (Node Incidences)','Fontsize',16,'Fontweight','Bold')
    legend('Topographic MaxEnt NM','NonTopographic MaxEnt NM','Topographic N&G NM','Nontopographic N&G NM','Weights Mat')
    xlabel('Row # of Unconstrained Weights Matrix','Fontsize',16,'Fontweight','Bold')
    ylabel('Sum of Row','Fontsize',16,'Fontweight','Bold')
    %
    figure, hold on, plot(sum(NMme_unc),'r'),plot(sum(NMng_unc),'m'),plot(sum(Wunc),'k')
    title('Nearest Neighbor Adjacency Row Sums (Node Incidences) Zoom [exclude Topographic Models]','Fontsize',16,'Fontweight','Bold')
    legend('Nontographic MaxEnt NM','Nontopographic N&G NM','Weights Mat')
    xlabel('Pixel Distance* (* distance in Weights Matrix is distance/orientation combination in image)','Fontsize',16,'Fontweight','Bold')
    ylabel('Expected Weight','Fontsize',16,'Fontweight','Bold')

    % (C). Plot node incidences (row sums) for MaxEnt for different values of max iteration with rowsums for weights matrix
    NMme20 = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,20,MEth,MEdyn,0);
    NMme50 = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,50,MEth,MEdyn,0);
    NMme100 = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,100,MEth,MEdyn,0);
    NMme200 = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,200,MEth,MEdyn,0);
    NMme500 = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,500,MEth,MEdyn,0);
    NMme1000 = compute_NM_MaxEnt(W,ones(size(W)),ximg,yimg,1000,MEth,MEdyn,0);
    figure, hold on, plot(sum(W),'k'), plot(sum(NMme20),'b'),plot(sum(NMme50),'r'),
    plot(sum(NMme100),'g'),plot(sum(NMme200),'m'),plot(sum(NMme500),'c'),plot(sum(NMme1000),'y'),
    title('Adjacency Row Sums (Node Incidences) vs. Topographic MaxEnt Row Sums','Fontsize',16,'Fontweight','Bold')
    legend('Adjacency (Truth)','TME maxIter 20','TME maxIter 50','TME maxIter 100','TME maxIter 200','TME maxIter 500','TME maxIter 1000')
    xlabel('Row #','Fontsize',16,'Fontweight','Bold')
    ylabel('Sum of Row','Fontsize',16,'Fontweight','Bold')




    %%
    % (A). Calculate pixel intensity depenedence on distance (and direction of
    % orientation because that is how W is built from rasterized image)
    [ B_W ] = PIV_vs_Dist(W,Wconn);
    [ B_tmeNM ] = PIV_vs_Dist(NMtme,Wconn);
    [ B_tngNM ] = PIV_vs_Dist(NMtng,Wconn);
    [ B_meNM ] = PIV_vs_Dist(NMme,ones(size(Wconn)));
    [ B_ngNM ] = PIV_vs_Dist(NMng,ones(size(Wconn)));

    % (B). Calculate pixel intensity depenedence on distance with unconstrained (no
    % explicit hardware constraints imposed by programmer)
    [ B_W_unc ] = PIV_vs_Dist(Wunc,ones(size(Wconn)));
    [ B_tmeNM_unc ] = PIV_vs_Dist(NMtme,ones(size(Wconn)));
    [ B_tngNM_unc ] = PIV_vs_Dist(NMtng,ones(size(Wconn)));
    [ B_meNM_unc ] = PIV_vs_Dist(NMme,ones(size(Wconn)));
    [ B_ngNM_unc ] = PIV_vs_Dist(NMng,ones(size(Wconn)));



    %%
    % (A). Plot Pixel Intensity Dependence on Distance (Diagonal Means) for Adjacency Matrix & different Null Models
    % Note: We can use different constraints (Nearest Neighbor, Radius of Interaction, None) to investigate image statistics.
    figure, hold on, 
    plot(B_tmeNM,'b','LineWidth',2),
    plot(B_meNM,'r','LineWidth',2),
    plot(B_tngNM,'g','LineWidth',2),
    plot(B_ngNM,'m','LineWidth',2),
    plot(B_W,'k','LineWidth',2)
    title('Nearest Neighbor Adjacency Diagonal Sums (Pixel Similarity vs. Distance/Orientation Relationship)','Fontsize',20,'Fontweight','Bold')
    legend('Topographic MaxEnt NM','NonTopographic MaxEnt NM','Topographic N&G NM','Nontopographic N&G NM','Weights Mat')
    xlabel('Pixel Distance* (* distance in Weights Matrix is distance/orientation combination in image)','Fontsize',16,'Fontweight','Bold')
    text(50,0.4,'<--- Topographical Models (blu & grn) and Weights (blk) right on top of one another.')
    ylabel('Expected Weight','Fontsize',16,'Fontweight','Bold')


    % (B). Plot Pixel Intensity Dependence on Distance (Diagonal Means) for Adjacency Matrix & different Null Models
    % Note: We can use different constraints (Nearest Neighbor, Radius of Interaction, None) to investigate image statistics.
    figure, hold on, 
    plot(B_tmeNM_unc,'b','LineWidth',2),
    plot(B_meNM_unc,'r','LineWidth',2),
    plot(B_tngNM_unc,'g','LineWidth',2),
    plot(B_ngNM_unc,'m','LineWidth',2),
    plot(B_W_unc,'k','LineWidth',2)
    title('Unconstrained Adjacency Diagonal Sums (Pixel Similarity vs. Distance/Orientation Relationship)','Fontsize',20,'Fontweight','Bold')
    legend('Topographic MaxEnt NM','NonTopographic MaxEnt NM','Topographic N&G NM','Nontopographic N&G NM','Weights Mat''Fontsize',16,'Fontweight','Bold')
    xlabel('Pixel Distance* (* distance in Weights Matrix is distance/orientation combination in image)','Fontsize',16,'Fontweight','Bold')
    ylabel('Expected Weight','Fontsize',16,'Fontweight','Bold')

    keyboard

end








