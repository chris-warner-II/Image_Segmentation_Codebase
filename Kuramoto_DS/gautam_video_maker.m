function gautam_video_maker(traceT,dims,ratio,im)
% this function takes in dynamical data (traceT), dimensions of array
% (dims), relative spacing in horizontal and vertical directions (ratio),
% and a background image (im), to produce a movie. This movie can be saved
% to a file.

subplot('position',[0 0 1 1]);
figure('Name','Traveling LFP waves in rat hippocampus');
if numel(dims) > 2
    probes = dims-min(dims(:))+1;%trace = trace(dims-min(dims(:)) + 1,:);
    dims = size(dims);
end
center = dims([2 1])/2;
scale = 20;
numPast = 500;
histT = nan*ones(numPast,2);
trange = std(traceT(:))*[-2 2];
% if exist('in','var') && ~isempty(in) && in
%     [xout yout] = meshgrid(ratio(2):ratio(2)*dims(2),ratio(1):ratio(1)*dims(1));
%     [xin yin] = meshgrid((1:dims(2))*ratio(2),(1:dims(1))*ratio(1));
%     dims = size(xout);
% end
for i = 1:size(traceT,2)%startInd + 500%
    temp = traceT(:,i);temp = temp(probes);
    %if exist('in','var') && ~isempty(in) && in
    %    temp = interp2(xin,yin,squeeze(temp),xout,yout,'cubic');
    %end
    temp = imfilter(temp,fspecial('gaussian',5,1));
    [xt yt] = angGradient(imresize(temp,1/2,'bilinear'));%myGradient(HT);%(reshape(HT(:,i),dims));
    [xs ys] = meshgrid(1:size(xt,2),1:size(xt,1));
    temp = real(temp);
    xt = -xt;yt = flipud(yt);tm = [mean(xt(:)) mean(yt(:))];
        histT = circshift(histT,[-1 0]);
        histT(numPast,:) = center([2 1]).*ratio+tm*scale;
        %subplot('Position',[0 0 1 .48]);
        %imagesc((1:dims(1))*ratio(1),(1:dims(2))*ratio(2),flipud(reshape(traceT(:,i),dims)),tRange*.8);hold on;axis off;
        %plot(histT(:,1),histT(:,2),'w','LineWidth',1.5);
        %plot(histT(:,1),histT(:,2),'k','LineWidth',1.5);
        %quiver(xt,yt,'k');
        %quiver(center(2)*ratio(1),center(1)*ratio(2),tm(1)*scale,tm(2)*scale,'w','LineWidth',5);
        %quiver(center(2)*ratio(1),center(1)*ratio(2),tm(1)*scale,tm(2)*scale,'k--','LineWidth',5);
        %hold off;
        %subplot('Position',[0 .52 1 .48]);%s1);%1,2,2);
        %temp = flipud(reshape(traceT(:,i),dims));
        colormap(jet);
        surf((1:dims(1))*ratio(1),(1:dims(2))*ratio(2),rot90(rot90(rot90(temp))),'edgecolor','none');set(gca,'zlim',trange);%hold on;
        title([num2str((i-1)*8/1250,2) ' s']);ylabel ('posterior-anterior');xlabel('ventral-dorsal');zlabel('voltage');
       % quiver3(center(2)*ratio(1),center(1)*ratio(2),trange(1),tm(1)*scale,tm(2)*scale,0,'k--','LineWidth',5);
       % plot3(histT(:,1),histT(:,2),trange(1)*ones(size(histT,1),1),'k','LineWidth',1.5);hold off;
        %quiver3(ys*ratio(1),xs*ratio(2),-1*ones(size(ys)),xt*scale,yt*scale,zeros(size(ys)),'k--','LineWidth',1);hold off;
        hold on;colormap(jet);
        surf([1 dims(1)]*ratio(1),[1 dims(2)]*ratio(2),repmat(-14, [2 2]),im,'facecolor','texture');
        set(gca,'xlim',[1 dims(1)]*ratio(1),'ylim',[1 dims(2)]*ratio(2),'zlim',[-14 14]);
        hold off;
        set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[],'fontsize',16);
%        imagesc((1:dims(1))*ratio(1),(1:dims(2))*ratio(2),flipud(reshape(trace(:,i),dims)),orRange*.8);axis off;
        axesLabelsAlign3D;drawnow;
        %m(i) = getframe(gcf);
end
%movie2avi(m,'crcnsSurf512.avi');

function fixStep(a,b)
global step
global h
step = round(get(h,'Value'));

