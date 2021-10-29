
% This script is to compare an oscillators actual (perturbed) frequency
% with its natural frequency to determine which oscillators it is
% influenced by / connected to.
%
% This is trying to compute a segmentation early by the initial transient
% dynamics of the oscillator network

cmap = colormap('jet');

osc = 1:12;


%% Look at Change in Frequency of Oscillator (actual freq compared to unperturbed/natural freq) 
h1=figure; hold on
h2=figure; hold on

for i = osc
    Yarg = diff(theta(:,i))./(tau*2*pi);
    Yarg = Yarg(Yarg>0); % = w(i);
    figure(h1), plot(Yarg,linestylee{i})
    figure(h2), plot(Yarg - w(i), linestylee{i})
    
    w_actual(i) = mean(Yarg);
end



%% Look at Pairwise Circular Distance (note this is same as below cosine distance)

for j = osc
    figure, hold on,
    k=0;
    for i = 1:12
        if(i~=j)
            k=k+1;
            plot(abs(circ_dist(theta(:,j),theta(:,i))),'Color',cmap(5*i,:),'LineWidth',2)
            xx(k) = mean( abs(circ_dist(theta(:,j),theta(:,i))) );
            yy(k) = std( abs(circ_dist(theta(:,j),theta(:,i))) );
            legos{k} = ['Osc ',num2str(i),' -- ',num2str(xx(k),2),' -- ',num2str(yy(k),2)];
        end
    end
    title(['Osc ',num2str(j),' Circular Distance'])
    legend(legos)
    xx
end



%% Pairwise Cosine Distance (note this is same as above circular distance)

% for j = 1:12
%     figure, hold on,
%     k=0;
%     for i = 1:12
%         if(i~=j)
%             k=k+1;
%             plot(abs(cos(theta(:,j) - theta(:,i))),'Color',cmap(5*i,:),'LineWidth',2)
%             xx(k) = mean( abs(cos(theta(:,j) - theta(:,i))) );
%             legos{k} = ['Osc ',num2str(i),' -- ',num2str(xx(k),2)];
%         end
%     end
%     title(['Osc ',num2str(j),' Circular Distance'])
%     legend(legos)
%     xx
% end