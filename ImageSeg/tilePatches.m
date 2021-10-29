function [patch,gTpatch] = tilePatches(image,groundTruth,pach,pnum)



patch = image(pach.ypbeg(pnum):pach.ypfin(pnum), pach.xpbeg(pnum):pach.xpfin(pnum));


for i = 1:numel(groundTruth)

    
    gTpatch{i} = groundTruth{i}.Segmentation(pach.ypbeg(pnum):pach.ypfin(pnum), pach.xpbeg(pnum):pach.xpfin(pnum));
    
    
end



patch = patch - min(patch(:));
patch = patch ./ max(patch(:)); % normalize patch to spread pixel values between 0 & 1.
