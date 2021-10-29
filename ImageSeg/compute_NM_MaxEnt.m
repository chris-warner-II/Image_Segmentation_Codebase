function [ NM, iter ] = compute_NM_MaxEnt(W,Wconn,ximg,yimg,maxIter,MEth,MEdyn,topo)

% syntax: [ NM ] = compute_null_model(W,Wconn,ximg,yimg,ds_fctr,maxIter,MEth,MEdyn)
%
% where NM is the output null model (expected connectivity of nodes in graph).
%       W is the weights (adjacency) matrix that describes actual connectivity in graph.
%      (ximg, yimg) are the dimensions of the actual image used to make W matrix.

r = 2; % number of bins to chop analog weights matrix into
n = size(W,1);
d = full(sum(W,2));  % degree sequence
if(~topo)
    Wconn = ones(size(Wconn));
    tag='Non-';
else
    tag='';
end
tic
[MEout,iter] = findMLE(d,r,zeros(n,1),maxIter,MEth,Wconn);  % find MLE
MEt = toc;

% in MaxEnt Max Likelihood Estimation, we init theta to all zeros (is this good)?

%% Visualize Dynamics of Max Entropy Algorithm with Movie
    if(MEdyn)
        disp(['Making movie of Max Ent Algorithm Dynamics'])
        tic
        iid = 1;
        Dyn = VideoWriter(['../output/MaxEnt/ME_',num2str(iid),'_Dyn_',num2str(ximg),...
            'x',num2str(yimg),'_imax',num2str(maxIter)],'Uncompressed AVI');
        Dyn.FrameRate = 3;
        open(Dyn); % open movie file
        for ii = 1:maxIter+1
            x1 = MEout(:,ii);
            hDyn = figure;
            imagesc(reshape(x1,ximg,yimg)/max(x1));
            colormap('bone')
            title([num2str(iid),' Max Ent Iter #',num2str(ii-1), ' : ',...
                num2str(round(MEt),'%d'), 'secs'],'FontSize',20,'FontWeight','Bold')
            axis off
            currFrame = getframe(hDyn);
            writeVideo(Dyn,currFrame);
            close(hDyn)
        end
        close(Dyn); % close movie file
        toc
        
        %imwrite(reshape(x1,ximg,yimg)/max(x1),'W232038.tif');
    end

Npot = MEout(:,size(MEout,2));

if(0)
figure, imagesc(reshape(Npot,yimg,ximg)),colormap('bone'),colorbar
title('MaxEnt Node Potentials','FontSize',20,'FontWeight','Bold')

end
NM = nodePots2edgeWaits(Npot);
NM=NM.*Wconn;

if ~topo
    xxx = mean(sum(NM) ./ sum(W)); % rescaling NM so that Rowsums match Rowsums of W matrix
    NM = NM./xxx;
end

if(0)
figure, imagesc(NM),colorbar, title(['MaxEnt ',tag,'Topographic Null Model'],'FontSize',20,'FontWeight','Bold')
end


%% Look into rescaling of unconverged MaxEnt Null Model values to better fit Adjacency (W) Row Sums
if(0)
MEout20 = findMLE(d,r,zeros(n,1),20,MEth,Wconn);  % find MLE
Npot20 = MEout20(:,size(MEout20,2));
NM20 = nodePots2edgeWaits(Npot20);
NM20=NM20.*Wconn;
xxx = mean(sum(NM20) ./ sum(W));
NM20 = NM20./xxx;


MEout500 = findMLE(d,r,zeros(n,1),500,MEth,Wconn);  % find MLE
Npot500 = MEout500(:,size(MEout500,2));
NM500 = nodePots2edgeWaits(Npot500);
NM500=NM500.*Wconn;

keyboard
end


%% Code to generate Full (Non-topological) MaxEnt Null Model to compare Topographical to Nontopo.
% MEout2 = findMLE(d,r,zeros(n,1),maxIter,MEth,ones(size(Wtopo)));  % find MLE
% Npot2 = MEout2(:,size(MEout2,2));
% figure, imagesc(reshape(Npot2,yimg,ximg)),colormap('bone'),colorbar
% title('MaxEnt Full Node Potentials','FontSize',20,'FontWeight','Bold')
% NM2 = nodePots2edgeWaits(Npot2);