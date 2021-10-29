function out = GaussianBox(x1,y1,x2,x3,y2,y3,mu_in,sig_in,mu_out,sig_out,disp)

% syntax: out = prob_dist_box(x1,y1,x2,x3,y2,y3,mu_in,sig_in,mu_out,sig_out,disp)
%
% This function makes an rectangular image of size (x1,y1).  The image
% contains a rectangle spanning from (x2-x3) and (y2-y3) within the
% rectangular image.  Both internal and external rectangles have pixel
% values picked from a gaussian probability distribution with mean Mu and
% standard deviation Sig.
%
% Would probably be useful to extend this to conditional probability
% distributions or 2D (x,y) joint distributions or something more.

%% define initial pixel values from probability distributions
rec_out = mu_out + sig_out.*randn(x1, y1);
rec_in = mu_in + sig_in.*randn(1+x3-x2, 1+y3-y2);


%% Iterate with while loop bc some may be < 0 or > 1 depending on mu & sigma
redo_out = 1;
while any(any(redo_out))
    redo_out = (rec_out < 0 | rec_out > 1);
    rec_out(redo_out) = 0;
    rec_out = rec_out + ( mu_out + sig_out.*randn(x1, y1) ) .* redo_out;
end
%
redo_in = 1;
while any(any(redo_in))
    redo_in = (rec_in < 0 | rec_in > 1);
    rec_in(redo_in) = 0;
    rec_in = rec_in + ( mu_in + sig_in.*randn(1+x3-x2, 1+y3-y2) ) .* redo_in;
end


%% put image of rectangle within rectangle together to be segmented
out = rec_out;
out(y2:y3,x2:x3) = rec_in;
out = out'; % flip because it tends to get flipped in making somehow.


%% Plot some images if disp flag is 1.
if disp
    % plot image to be segmented
    figure, imagesc(out), colormap('bone'), colorbar

    % Plot histograms of inside and outside rectangle pixel value distributions
    [h_out,x_out] = hist(rec_out(:),1000);
    [h_in,x_in] = hist(rec_in(:),1000);
    figure, hold on
    plot(x_out,h_out/max(h_out),'r'), plot(x_in,h_in/max(h_in),'b');
end
