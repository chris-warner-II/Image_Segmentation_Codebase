
%% To run image_seg on 5 image patches for {IsoDiff, AA, Mod_SKH} & rM={1,3,10}
if(0)
    
    Loop_ImgSegMethodsD(12,[1,3,10],0.2,{[51,51]},'BSDS_patch',[24,24],0);
    Loop_ImgSegMethodsD(3,[1,3,10],0.2,{[51,51]},'BSDS_patch',[24,24],0);
    Loop_ImgSegMethodsD(7,[1,3,10],0.2,{[51,51]},'BSDS_patch',[24,24],0);
    %
    Loop_ImgSegMethodsD(12,[1,3,10],0.2,{[51,51]},'BSDS_patch',[37,37],0);
    Loop_ImgSegMethodsD(3,[1,3,10],0.2,{[51,51]},'BSDS_patch',[37,37],0);
    Loop_ImgSegMethodsD(7,[1,3,10],0.2,{[51,51]},'BSDS_patch',[37,37],0);
    %
    Loop_ImgSegMethodsD(12,[1,3,10],0.2,{[51,51]},'BSDS_patch',[46,46],0);
    Loop_ImgSegMethodsD(3,[1,3,10],0.2,{[51,51]},'BSDS_patch',[46,46],0);
    Loop_ImgSegMethodsD(7,[1,3,10],0.2,{[51,51]},'BSDS_patch',[46,46],0);
    %
    Loop_ImgSegMethodsD(12,[1,3,10],0.2,{[51,51]},'BSDS_patch',[54,54],0);
    Loop_ImgSegMethodsD(3,[1,3,10],0.2,{[51,51]},'BSDS_patch',[54,54],0);
    Loop_ImgSegMethodsD(7,[1,3,10],0.2,{[51,51]},'BSDS_patch',[54,54],0);
    
end

%% To run plot_phaseAtClk_vs_Kscale for these other 5 images.
if(1)
    
    img_ptch = {'100080_ptch2','10081_ptch1','101027_ptch3','102062_ptch1','103041_ptch1','103078_ptch3'};

    for i = 1:numel(img_ptch)
        plot_phaseAtClk_vs_Kscale('IsoDiff', img_ptch{i}, 'rM1_');
        plot_phaseAtClk_vs_Kscale('IsoDiff', img_ptch{i}, 'rM3_');
        plot_phaseAtClk_vs_Kscale('IsoDiff', img_ptch{i}, 'rM10_');
        %
        plot_phaseAtClk_vs_Kscale('AAnrm', img_ptch{i}, 'rM1_');
        plot_phaseAtClk_vs_Kscale('AAnrm', img_ptch{i}, 'rM3_');
        plot_phaseAtClk_vs_Kscale('AAnrm', img_ptch{i}, 'rM10_');
        %
        plot_phaseAtClk_vs_Kscale('Mod_SKHAdj', img_ptch{i}, 'rM1_');
        plot_phaseAtClk_vs_Kscale('Mod_SKHAdj', img_ptch{i}, 'rM3_');
        plot_phaseAtClk_vs_Kscale('Mod_SKHAdj', img_ptch{i}, 'rM10_');
    end

end


%% Make noise to alert me its done.
while(1)
    beep
    pause(1)
end