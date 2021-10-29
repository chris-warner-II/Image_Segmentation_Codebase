function compare_methods_AUC_ROC_resC

    % syntax:  compare_methods_AUC_ROC_resB(fileType, fileSize)
    %
    %  This function will spit out plots of average Area under ROC curve
    %  with STD error bars for the 4 different netMethods and for 4
    %  different segMethods also.
    %
    % fileType = 'BSDS_patch';
    % fileSize = '51x51_ds1';




    %% Load in mat file and set up directories to different data.
    [dirPre,sizeGoodIm] = onCluster;

%     matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/AUC_ROC_results/'];
%     
%     imgFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/AUC_ROC_results/'];
%     
%     do_AA = 1;
%     do_GL = 1;
%     do_SK = 1;
%     do_NG = 1;
%     do_Iso = 1;
    



    if(0)

    
        %% Load Results from smaller (51x51) Image Patches

        AA_51_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_51_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_51_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        AA_51_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_51_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_51_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        AA_51_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_51_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_51_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        AA_51_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_51_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_51_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        GL_51_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_51_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_51_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        GL_51_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_51_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_51_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        GL_51_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_51_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_51_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        GL_51_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_51_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_51_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        SK_51_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_51_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_51_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        SK_51_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_51_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_51_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        SK_51_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_51_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_51_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        SK_51_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_51_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_51_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        NG_51_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_51_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_51_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        NG_51_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_51_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_51_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        NG_51_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_51_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_51_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        NG_51_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_51_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_51_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        Iso_51_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM1_NF_60_0_kssml.mat']);
        Iso_51_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM1_NF_60_0_ksmid.mat']);
        Iso_51_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM1_NF_60_0_kslrg.mat']);
        %
        Iso_51_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM3_NF_60_0_kssml.mat']);
        Iso_51_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM3_NF_60_0_ksmid.mat']);
        Iso_51_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM3_NF_60_0_kslrg.mat']);
        %
        Iso_51_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM5_NF_60_0_kssml.mat']);
        Iso_51_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM5_NF_60_0_ksmid.mat']);
        Iso_51_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM5_NF_60_0_kslrg.mat']);
        %
        Iso_51_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM10_NF_60_0_kssml.mat']);
        Iso_51_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM10_NF_60_0_ksmid.mat']);
        Iso_51_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM10_NF_60_0_kslrg.mat']);
        %

        %% Load Results from larger (101x101) Image Patches

        AA_101_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_101_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_101_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        AA_101_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_101_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_101_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        AA_101_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_101_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_101_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        AA_101_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        AA_101_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        AA_101_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/AAnrm/AUCdata_AAnrm_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        GL_101_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_101_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_101_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        GL_101_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_101_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_101_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        GL_101_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_101_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_101_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        GL_101_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        GL_101_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        GL_101_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/GLnrm/AUCdata_GLnrm_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        SK_101_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_101_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_101_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        SK_101_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_101_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_101_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        SK_101_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_101_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_101_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        SK_101_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        SK_101_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        SK_101_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        NG_101_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM1_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_101_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM1_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_101_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        NG_101_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM3_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_101_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_101_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        NG_101_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM5_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_101_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM5_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_101_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %
        NG_101_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM10_sDInf_sP0p2_NF_60_0_kssml.mat']);
        NG_101_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM10_sDInf_sP0p2_NF_60_0_ksmid.mat']);
        NG_101_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/Mod_N&G/AUCdata_Mod_N&G_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %

        Iso_101_1_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM1_NF_60_0_kssml.mat']);
        Iso_101_1_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM1_NF_60_0_ksmid.mat']);
        Iso_101_1_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM1_NF_60_0_kslrg.mat']);
        %
        Iso_101_3_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM3_NF_60_0_kssml.mat']);
        Iso_101_3_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM3_NF_60_0_ksmid.mat']);
        Iso_101_3_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM3_NF_60_0_kslrg.mat']);
        %
        Iso_101_5_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM5_NF_60_0_kssml.mat']);
        Iso_101_5_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM5_NF_60_0_ksmid.mat']);
        Iso_101_5_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM5_NF_60_0_kslrg.mat']);
        %
        Iso_101_10_sml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM10_NF_60_0_kssml.mat']);
        Iso_101_10_mid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM10_NF_60_0_ksmid.mat']);
        Iso_101_10_lrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/AUC_ROC_results/IsoDiff/AUCdata_IsoDiff_rM10_NF_60_0_kslrg.mat']);
        %






        %% Histogram, Mean & Std for Kur Upper Bound for Isotropic Diffusion Methods

        bins = [-0.5:0.05:0.5];

        % histograms for smaller (51x51) image patches
        hIso_51_1_lrg = hist(Iso_51_1_lrg.kur_ROC_AUC_GS_accum - Iso_51_1_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_51_1_mid = hist(Iso_51_1_mid.kur_ROC_AUC_GS_accum - Iso_51_1_mid.im_ROC_AUC_1D_accum, bins);
        hIso_51_1_sml = hist(Iso_51_1_sml.kur_ROC_AUC_GS_accum - Iso_51_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hIso_51_3_lrg = hist(Iso_51_3_lrg.kur_ROC_AUC_GS_accum - Iso_51_3_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_51_3_mid = hist(Iso_51_3_mid.kur_ROC_AUC_GS_accum - Iso_51_3_mid.im_ROC_AUC_1D_accum, bins);
        hIso_51_3_sml = hist(Iso_51_3_sml.kur_ROC_AUC_GS_accum - Iso_51_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hIso_51_5_lrg = hist(Iso_51_5_lrg.kur_ROC_AUC_GS_accum - Iso_51_5_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_51_5_mid = hist(Iso_51_5_mid.kur_ROC_AUC_GS_accum - Iso_51_5_mid.im_ROC_AUC_1D_accum, bins);
        hIso_51_5_sml = hist(Iso_51_5_sml.kur_ROC_AUC_GS_accum - Iso_51_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hIso_51_10_lrg = hist(Iso_51_10_lrg.kur_ROC_AUC_GS_accum - Iso_51_10_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_51_10_mid = hist(Iso_51_10_mid.kur_ROC_AUC_GS_accum - Iso_51_10_mid.im_ROC_AUC_1D_accum, bins);
        hIso_51_10_sml = hist(Iso_51_10_sml.kur_ROC_AUC_GS_accum - Iso_51_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (51x51) image patches
        mIso_51_1_lrg = mean(Iso_51_1_lrg.kur_ROC_AUC_GS_accum - Iso_51_1_lrg.im_ROC_AUC_1D_accum);
        mIso_51_1_mid = mean(Iso_51_1_mid.kur_ROC_AUC_GS_accum - Iso_51_1_mid.im_ROC_AUC_1D_accum);
        mIso_51_1_sml = mean(Iso_51_1_sml.kur_ROC_AUC_GS_accum - Iso_51_1_sml.im_ROC_AUC_1D_accum);
        %
        mIso_51_3_lrg = mean(Iso_51_3_lrg.kur_ROC_AUC_GS_accum - Iso_51_3_lrg.im_ROC_AUC_1D_accum);
        mIso_51_3_mid = mean(Iso_51_3_mid.kur_ROC_AUC_GS_accum - Iso_51_3_mid.im_ROC_AUC_1D_accum);
        mIso_51_3_sml = mean(Iso_51_3_sml.kur_ROC_AUC_GS_accum - Iso_51_3_sml.im_ROC_AUC_1D_accum);
        %
        mIso_51_5_lrg = mean(Iso_51_5_lrg.kur_ROC_AUC_GS_accum - Iso_51_5_lrg.im_ROC_AUC_1D_accum);
        mIso_51_5_mid = mean(Iso_51_5_mid.kur_ROC_AUC_GS_accum - Iso_51_5_mid.im_ROC_AUC_1D_accum);
        mIso_51_5_sml = mean(Iso_51_5_sml.kur_ROC_AUC_GS_accum - Iso_51_5_sml.im_ROC_AUC_1D_accum);
        %
        mIso_51_10_lrg = mean(Iso_51_10_lrg.kur_ROC_AUC_GS_accum - Iso_51_10_lrg.im_ROC_AUC_1D_accum);
        mIso_51_10_mid = mean(Iso_51_10_mid.kur_ROC_AUC_GS_accum - Iso_51_10_mid.im_ROC_AUC_1D_accum);
        mIso_51_10_sml = mean(Iso_51_10_sml.kur_ROC_AUC_GS_accum - Iso_51_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (51x51) image patches
        sIso_51_1_lrg = std(Iso_51_1_lrg.kur_ROC_AUC_GS_accum - Iso_51_1_lrg.im_ROC_AUC_1D_accum);
        sIso_51_1_mid = std(Iso_51_1_mid.kur_ROC_AUC_GS_accum - Iso_51_1_mid.im_ROC_AUC_1D_accum);
        sIso_51_1_sml = std(Iso_51_1_sml.kur_ROC_AUC_GS_accum - Iso_51_1_sml.im_ROC_AUC_1D_accum);
        %
        sIso_51_3_lrg = std(Iso_51_3_lrg.kur_ROC_AUC_GS_accum - Iso_51_3_lrg.im_ROC_AUC_1D_accum);
        sIso_51_3_mid = std(Iso_51_3_mid.kur_ROC_AUC_GS_accum - Iso_51_3_mid.im_ROC_AUC_1D_accum);
        sIso_51_3_sml = std(Iso_51_3_sml.kur_ROC_AUC_GS_accum - Iso_51_3_sml.im_ROC_AUC_1D_accum);
        %
        sIso_51_5_lrg = std(Iso_51_5_lrg.kur_ROC_AUC_GS_accum - Iso_51_5_lrg.im_ROC_AUC_1D_accum);
        sIso_51_5_mid = std(Iso_51_5_mid.kur_ROC_AUC_GS_accum - Iso_51_5_mid.im_ROC_AUC_1D_accum);
        sIso_51_5_sml = std(Iso_51_5_sml.kur_ROC_AUC_GS_accum - Iso_51_5_sml.im_ROC_AUC_1D_accum);
        %
        sIso_51_10_lrg = std(Iso_51_10_lrg.kur_ROC_AUC_GS_accum - Iso_51_10_lrg.im_ROC_AUC_1D_accum);
        sIso_51_10_mid = std(Iso_51_10_mid.kur_ROC_AUC_GS_accum - Iso_51_10_mid.im_ROC_AUC_1D_accum);
        sIso_51_10_sml = std(Iso_51_10_sml.kur_ROC_AUC_GS_accum - Iso_51_10_sml.im_ROC_AUC_1D_accum);


        % histograms for larger (101x101) image patches
        hIso_101_1_lrg = hist(Iso_101_1_lrg.kur_ROC_AUC_GS_accum - Iso_101_1_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_101_1_mid = hist(Iso_101_1_mid.kur_ROC_AUC_GS_accum - Iso_101_1_mid.im_ROC_AUC_1D_accum, bins);
        hIso_101_1_sml = hist(Iso_101_1_sml.kur_ROC_AUC_GS_accum - Iso_101_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hIso_101_3_lrg = hist(Iso_101_3_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_101_3_mid = hist(Iso_101_3_mid.kur_ROC_AUC_GS_accum - Iso_101_3_mid.im_ROC_AUC_1D_accum, bins);
        hIso_101_3_sml = hist(Iso_101_3_sml.kur_ROC_AUC_GS_accum - Iso_101_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hIso_101_5_lrg = hist(Iso_101_5_lrg.kur_ROC_AUC_GS_accum - Iso_101_5_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_101_5_mid = hist(Iso_101_5_mid.kur_ROC_AUC_GS_accum - Iso_101_5_mid.im_ROC_AUC_1D_accum, bins);
        hIso_101_5_sml = hist(Iso_101_5_sml.kur_ROC_AUC_GS_accum - Iso_101_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hIso_101_10_lrg = hist(Iso_101_10_lrg.kur_ROC_AUC_GS_accum - Iso_101_10_lrg.im_ROC_AUC_1D_accum, bins);
        hIso_101_10_mid = hist(Iso_101_10_mid.kur_ROC_AUC_GS_accum - Iso_101_10_mid.im_ROC_AUC_1D_accum, bins);
        hIso_101_10_sml = hist(Iso_101_10_sml.kur_ROC_AUC_GS_accum - Iso_101_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (101x101) image patches
        mIso_101_1_lrg = mean(Iso_101_1_lrg.kur_ROC_AUC_GS_accum - Iso_101_1_lrg.im_ROC_AUC_1D_accum);
        mIso_101_1_mid = mean(Iso_101_1_mid.kur_ROC_AUC_GS_accum - Iso_101_1_mid.im_ROC_AUC_1D_accum);
        mIso_101_1_sml = mean(Iso_101_1_sml.kur_ROC_AUC_GS_accum - Iso_101_1_sml.im_ROC_AUC_1D_accum);
        %
        mIso_101_3_lrg = mean(Iso_101_3_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_lrg.im_ROC_AUC_1D_accum);
        mIso_101_3_mid = mean(Iso_101_3_mid.kur_ROC_AUC_GS_accum - Iso_101_3_mid.im_ROC_AUC_1D_accum);
        mIso_101_3_sml = mean(Iso_101_3_sml.kur_ROC_AUC_GS_accum - Iso_101_3_sml.im_ROC_AUC_1D_accum);
        %
        mIso_101_5_lrg = mean(Iso_101_5_lrg.kur_ROC_AUC_GS_accum - Iso_101_5_lrg.im_ROC_AUC_1D_accum);
        mIso_101_5_mid = mean(Iso_101_5_mid.kur_ROC_AUC_GS_accum - Iso_101_5_mid.im_ROC_AUC_1D_accum);
        mIso_101_5_sml = mean(Iso_101_5_sml.kur_ROC_AUC_GS_accum - Iso_101_5_sml.im_ROC_AUC_1D_accum);
        %
        mIso_101_10_lrg = mean(Iso_101_10_lrg.kur_ROC_AUC_GS_accum - Iso_101_10_lrg.im_ROC_AUC_1D_accum);
        mIso_101_10_mid = mean(Iso_101_10_mid.kur_ROC_AUC_GS_accum - Iso_101_10_mid.im_ROC_AUC_1D_accum);
        mIso_101_10_sml = mean(Iso_101_10_sml.kur_ROC_AUC_GS_accum - Iso_101_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (101x101) image patches
        sIso_101_1_lrg = std(Iso_101_1_lrg.kur_ROC_AUC_GS_accum - Iso_101_1_lrg.im_ROC_AUC_1D_accum);
        sIso_101_1_mid = std(Iso_101_1_mid.kur_ROC_AUC_GS_accum - Iso_101_1_mid.im_ROC_AUC_1D_accum);
        sIso_101_1_sml = std(Iso_101_1_sml.kur_ROC_AUC_GS_accum - Iso_101_1_sml.im_ROC_AUC_1D_accum);
        %
        sIso_101_3_lrg = std(Iso_101_3_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_lrg.im_ROC_AUC_1D_accum);
        sIso_101_3_mid = std(Iso_101_3_mid.kur_ROC_AUC_GS_accum - Iso_101_3_mid.im_ROC_AUC_1D_accum);
        sIso_101_3_sml = std(Iso_101_3_sml.kur_ROC_AUC_GS_accum - Iso_101_3_sml.im_ROC_AUC_1D_accum);
        %
        sIso_101_5_lrg = std(Iso_101_5_lrg.kur_ROC_AUC_GS_accum - Iso_101_5_lrg.im_ROC_AUC_1D_accum);
        sIso_101_5_mid = std(Iso_101_5_mid.kur_ROC_AUC_GS_accum - Iso_101_5_mid.im_ROC_AUC_1D_accum);
        sIso_101_5_sml = std(Iso_101_5_sml.kur_ROC_AUC_GS_accum - Iso_101_5_sml.im_ROC_AUC_1D_accum);
        %
        sIso_101_10_lrg = std(Iso_101_10_lrg.kur_ROC_AUC_GS_accum - Iso_101_10_lrg.im_ROC_AUC_1D_accum);
        sIso_101_10_mid = std(Iso_101_10_mid.kur_ROC_AUC_GS_accum - Iso_101_10_mid.im_ROC_AUC_1D_accum);
        sIso_101_10_sml = std(Iso_101_10_sml.kur_ROC_AUC_GS_accum - Iso_101_10_sml.im_ROC_AUC_1D_accum);





        %% Histogram, Mean & Std for Kur Upper Bound for Graph Laplacian Methods

        % histograms for smaller (51x51) image patches
        hGL_51_1_lrg = hist(GL_51_1_lrg.kur_ROC_AUC_GS_accum - GL_51_1_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_51_1_mid = hist(GL_51_1_mid.kur_ROC_AUC_GS_accum - GL_51_1_mid.im_ROC_AUC_1D_accum, bins);
        hGL_51_1_sml = hist(GL_51_1_sml.kur_ROC_AUC_GS_accum - GL_51_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGL_51_3_lrg = hist(GL_51_3_lrg.kur_ROC_AUC_GS_accum - GL_51_3_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_51_3_mid = hist(GL_51_3_mid.kur_ROC_AUC_GS_accum - GL_51_3_mid.im_ROC_AUC_1D_accum, bins);
        hGL_51_3_sml = hist(GL_51_3_sml.kur_ROC_AUC_GS_accum - GL_51_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGL_51_5_lrg = hist(GL_51_5_lrg.kur_ROC_AUC_GS_accum - GL_51_5_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_51_5_mid = hist(GL_51_5_mid.kur_ROC_AUC_GS_accum - GL_51_5_mid.im_ROC_AUC_1D_accum, bins);
        hGL_51_5_sml = hist(GL_51_5_sml.kur_ROC_AUC_GS_accum - GL_51_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGL_51_10_lrg = hist(GL_51_10_lrg.kur_ROC_AUC_GS_accum - GL_51_10_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_51_10_mid = hist(GL_51_10_mid.kur_ROC_AUC_GS_accum - GL_51_10_mid.im_ROC_AUC_1D_accum, bins);
        hGL_51_10_sml = hist(GL_51_10_sml.kur_ROC_AUC_GS_accum - GL_51_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (51x51) image patches
        mGL_51_1_lrg = mean(GL_51_1_lrg.kur_ROC_AUC_GS_accum - GL_51_1_lrg.im_ROC_AUC_1D_accum);
        mGL_51_1_mid = mean(GL_51_1_mid.kur_ROC_AUC_GS_accum - GL_51_1_mid.im_ROC_AUC_1D_accum);
        mGL_51_1_sml = mean(GL_51_1_sml.kur_ROC_AUC_GS_accum - GL_51_1_sml.im_ROC_AUC_1D_accum);
        %
        mGL_51_3_lrg = mean(GL_51_3_lrg.kur_ROC_AUC_GS_accum - GL_51_3_lrg.im_ROC_AUC_1D_accum);
        mGL_51_3_mid = mean(GL_51_3_mid.kur_ROC_AUC_GS_accum - GL_51_3_mid.im_ROC_AUC_1D_accum);
        mGL_51_3_sml = mean(GL_51_3_sml.kur_ROC_AUC_GS_accum - GL_51_3_sml.im_ROC_AUC_1D_accum);
        %
        mGL_51_5_lrg = mean(GL_51_5_lrg.kur_ROC_AUC_GS_accum - GL_51_5_lrg.im_ROC_AUC_1D_accum);
        mGL_51_5_mid = mean(GL_51_5_mid.kur_ROC_AUC_GS_accum - GL_51_5_mid.im_ROC_AUC_1D_accum);
        mGL_51_5_sml = mean(GL_51_5_sml.kur_ROC_AUC_GS_accum - GL_51_5_sml.im_ROC_AUC_1D_accum);
        %
        mGL_51_10_lrg = mean(GL_51_10_lrg.kur_ROC_AUC_GS_accum - GL_51_10_lrg.im_ROC_AUC_1D_accum);
        mGL_51_10_mid = mean(GL_51_10_mid.kur_ROC_AUC_GS_accum - GL_51_10_mid.im_ROC_AUC_1D_accum);
        mGL_51_10_sml = mean(GL_51_10_sml.kur_ROC_AUC_GS_accum - GL_51_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (51x51) image patches
        sGL_51_1_lrg = std(GL_51_1_lrg.kur_ROC_AUC_GS_accum - GL_51_1_lrg.im_ROC_AUC_1D_accum);
        sGL_51_1_mid = std(GL_51_1_mid.kur_ROC_AUC_GS_accum - GL_51_1_mid.im_ROC_AUC_1D_accum);
        sGL_51_1_sml = std(GL_51_1_sml.kur_ROC_AUC_GS_accum - GL_51_1_sml.im_ROC_AUC_1D_accum);
        %
        sGL_51_3_lrg = std(GL_51_3_lrg.kur_ROC_AUC_GS_accum - GL_51_3_lrg.im_ROC_AUC_1D_accum);
        sGL_51_3_mid = std(GL_51_3_mid.kur_ROC_AUC_GS_accum - GL_51_3_mid.im_ROC_AUC_1D_accum);
        sGL_51_3_sml = std(GL_51_3_sml.kur_ROC_AUC_GS_accum - GL_51_3_sml.im_ROC_AUC_1D_accum);
        %
        sGL_51_5_lrg = std(GL_51_5_lrg.kur_ROC_AUC_GS_accum - GL_51_5_lrg.im_ROC_AUC_1D_accum);
        sGL_51_5_mid = std(GL_51_5_mid.kur_ROC_AUC_GS_accum - GL_51_5_mid.im_ROC_AUC_1D_accum);
        sGL_51_5_sml = std(GL_51_5_sml.kur_ROC_AUC_GS_accum - GL_51_5_sml.im_ROC_AUC_1D_accum);
        %
        sGL_51_10_lrg = std(GL_51_10_lrg.kur_ROC_AUC_GS_accum - GL_51_10_lrg.im_ROC_AUC_1D_accum);
        sGL_51_10_mid = std(GL_51_10_mid.kur_ROC_AUC_GS_accum - GL_51_10_mid.im_ROC_AUC_1D_accum);
        sGL_51_10_sml = std(GL_51_10_sml.kur_ROC_AUC_GS_accum - GL_51_10_sml.im_ROC_AUC_1D_accum);


        % histograms for larger (101x101) image patches
        hGL_101_1_lrg = hist(GL_101_1_lrg.kur_ROC_AUC_GS_accum - GL_101_1_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_101_1_mid = hist(GL_101_1_mid.kur_ROC_AUC_GS_accum - GL_101_1_mid.im_ROC_AUC_1D_accum, bins);
        hGL_101_1_sml = hist(GL_101_1_sml.kur_ROC_AUC_GS_accum - GL_101_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGL_101_3_lrg = hist(GL_101_3_lrg.kur_ROC_AUC_GS_accum - GL_101_3_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_101_3_mid = hist(GL_101_3_mid.kur_ROC_AUC_GS_accum - GL_101_3_mid.im_ROC_AUC_1D_accum, bins);
        hGL_101_3_sml = hist(GL_101_3_sml.kur_ROC_AUC_GS_accum - GL_101_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGL_101_5_lrg = hist(GL_101_5_lrg.kur_ROC_AUC_GS_accum - GL_101_5_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_101_5_mid = hist(GL_101_5_mid.kur_ROC_AUC_GS_accum - GL_101_5_mid.im_ROC_AUC_1D_accum, bins);
        hGL_101_5_sml = hist(GL_101_5_sml.kur_ROC_AUC_GS_accum - GL_101_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGL_101_10_lrg = hist(GL_101_10_lrg.kur_ROC_AUC_GS_accum - GL_101_10_lrg.im_ROC_AUC_1D_accum, bins);
        hGL_101_10_mid = hist(GL_101_10_mid.kur_ROC_AUC_GS_accum - GL_101_10_mid.im_ROC_AUC_1D_accum, bins);
        hGL_101_10_sml = hist(GL_101_10_sml.kur_ROC_AUC_GS_accum - GL_101_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (101x101) image patches
        mGL_101_1_lrg = mean(GL_101_1_lrg.kur_ROC_AUC_GS_accum - GL_101_1_lrg.im_ROC_AUC_1D_accum);
        mGL_101_1_mid = mean(GL_101_1_mid.kur_ROC_AUC_GS_accum - GL_101_1_mid.im_ROC_AUC_1D_accum);
        mGL_101_1_sml = mean(GL_101_1_sml.kur_ROC_AUC_GS_accum - GL_101_1_sml.im_ROC_AUC_1D_accum);
        %
        mGL_101_3_lrg = mean(GL_101_3_lrg.kur_ROC_AUC_GS_accum - GL_101_3_lrg.im_ROC_AUC_1D_accum);
        mGL_101_3_mid = mean(GL_101_3_mid.kur_ROC_AUC_GS_accum - GL_101_3_mid.im_ROC_AUC_1D_accum);
        mGL_101_3_sml = mean(GL_101_3_sml.kur_ROC_AUC_GS_accum - GL_101_3_sml.im_ROC_AUC_1D_accum);
        %
        mGL_101_5_lrg = mean(GL_101_5_lrg.kur_ROC_AUC_GS_accum - GL_101_5_lrg.im_ROC_AUC_1D_accum);
        mGL_101_5_mid = mean(GL_101_5_mid.kur_ROC_AUC_GS_accum - GL_101_5_mid.im_ROC_AUC_1D_accum);
        mGL_101_5_sml = mean(GL_101_5_sml.kur_ROC_AUC_GS_accum - GL_101_5_sml.im_ROC_AUC_1D_accum);
        %
        mGL_101_10_lrg = mean(GL_101_10_lrg.kur_ROC_AUC_GS_accum - GL_101_10_lrg.im_ROC_AUC_1D_accum);
        mGL_101_10_mid = mean(GL_101_10_mid.kur_ROC_AUC_GS_accum - GL_101_10_mid.im_ROC_AUC_1D_accum);
        mGL_101_10_sml = mean(GL_101_10_sml.kur_ROC_AUC_GS_accum - GL_101_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (101x101) image patches
        sGL_101_1_lrg = std(GL_101_1_lrg.kur_ROC_AUC_GS_accum - GL_101_1_lrg.im_ROC_AUC_1D_accum);
        sGL_101_1_mid = std(GL_101_1_mid.kur_ROC_AUC_GS_accum - GL_101_1_mid.im_ROC_AUC_1D_accum);
        sGL_101_1_sml = std(GL_101_1_sml.kur_ROC_AUC_GS_accum - GL_101_1_sml.im_ROC_AUC_1D_accum);
        %
        sGL_101_3_lrg = std(GL_101_3_lrg.kur_ROC_AUC_GS_accum - GL_101_3_lrg.im_ROC_AUC_1D_accum);
        sGL_101_3_mid = std(GL_101_3_mid.kur_ROC_AUC_GS_accum - GL_101_3_mid.im_ROC_AUC_1D_accum);
        sGL_101_3_sml = std(GL_101_3_sml.kur_ROC_AUC_GS_accum - GL_101_3_sml.im_ROC_AUC_1D_accum);
        %
        sGL_101_5_lrg = std(GL_101_5_lrg.kur_ROC_AUC_GS_accum - GL_101_5_lrg.im_ROC_AUC_1D_accum);
        sGL_101_5_mid = std(GL_101_5_mid.kur_ROC_AUC_GS_accum - GL_101_5_mid.im_ROC_AUC_1D_accum);
        sGL_101_5_sml = std(GL_101_5_sml.kur_ROC_AUC_GS_accum - GL_101_5_sml.im_ROC_AUC_1D_accum);
        %
        sGL_101_10_lrg = std(GL_101_10_lrg.kur_ROC_AUC_GS_accum - GL_101_10_lrg.im_ROC_AUC_1D_accum);
        sGL_101_10_mid = std(GL_101_10_mid.kur_ROC_AUC_GS_accum - GL_101_10_mid.im_ROC_AUC_1D_accum);
        sGL_101_10_sml = std(GL_101_10_sml.kur_ROC_AUC_GS_accum - GL_101_10_sml.im_ROC_AUC_1D_accum);



        %% Histogram, Mean & Std for Kur Upper Bound for Graph Laplacian Methods

        % histograms for smaller (51x51) image patches
        hGLev_51_1_lrg = hist(GL_51_1_lrg.ev2o_ROC_AUC_1D_accum - GL_51_1_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_51_1_mid = hist(GL_51_1_mid.ev2o_ROC_AUC_1D_accum - GL_51_1_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_51_1_sml = hist(GL_51_1_sml.ev2o_ROC_AUC_1D_accum - GL_51_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGLev_51_3_lrg = hist(GL_51_3_lrg.ev2o_ROC_AUC_1D_accum - GL_51_3_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_51_3_mid = hist(GL_51_3_mid.ev2o_ROC_AUC_1D_accum - GL_51_3_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_51_3_sml = hist(GL_51_3_sml.ev2o_ROC_AUC_1D_accum - GL_51_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGLev_51_5_lrg = hist(GL_51_5_lrg.ev2o_ROC_AUC_1D_accum - GL_51_5_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_51_5_mid = hist(GL_51_5_mid.ev2o_ROC_AUC_1D_accum - GL_51_5_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_51_5_sml = hist(GL_51_5_sml.ev2o_ROC_AUC_1D_accum - GL_51_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGLev_51_10_lrg = hist(GL_51_10_lrg.ev2o_ROC_AUC_1D_accum - GL_51_10_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_51_10_mid = hist(GL_51_10_mid.ev2o_ROC_AUC_1D_accum - GL_51_10_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_51_10_sml = hist(GL_51_10_sml.ev2o_ROC_AUC_1D_accum - GL_51_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (51x51) image patches
        mGLev_51_1_lrg = mean(GL_51_1_lrg.ev2o_ROC_AUC_1D_accum - GL_51_1_lrg.im_ROC_AUC_1D_accum);
        mGLev_51_1_mid = mean(GL_51_1_mid.ev2o_ROC_AUC_1D_accum - GL_51_1_mid.im_ROC_AUC_1D_accum);
        mGLev_51_1_sml = mean(GL_51_1_sml.ev2o_ROC_AUC_1D_accum - GL_51_1_sml.im_ROC_AUC_1D_accum);
        %
        mGLev_51_3_lrg = mean(GL_51_3_lrg.ev2o_ROC_AUC_1D_accum - GL_51_3_lrg.im_ROC_AUC_1D_accum);
        mGLev_51_3_mid = mean(GL_51_3_mid.ev2o_ROC_AUC_1D_accum - GL_51_3_mid.im_ROC_AUC_1D_accum);
        mGLev_51_3_sml = mean(GL_51_3_sml.ev2o_ROC_AUC_1D_accum - GL_51_3_sml.im_ROC_AUC_1D_accum);
        %
        mGLev_51_5_lrg = mean(GL_51_5_lrg.ev2o_ROC_AUC_1D_accum - GL_51_5_lrg.im_ROC_AUC_1D_accum);
        mGLev_51_5_mid = mean(GL_51_5_mid.ev2o_ROC_AUC_1D_accum - GL_51_5_mid.im_ROC_AUC_1D_accum);
        mGLev_51_5_sml = mean(GL_51_5_sml.ev2o_ROC_AUC_1D_accum - GL_51_5_sml.im_ROC_AUC_1D_accum);
        %
        mGLev_51_10_lrg = mean(GL_51_10_lrg.ev2o_ROC_AUC_1D_accum - GL_51_10_lrg.im_ROC_AUC_1D_accum);
        mGLev_51_10_mid = mean(GL_51_10_mid.ev2o_ROC_AUC_1D_accum - GL_51_10_mid.im_ROC_AUC_1D_accum);
        mGLev_51_10_sml = mean(GL_51_10_sml.ev2o_ROC_AUC_1D_accum - GL_51_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (51x51) image patches
        sGLev_51_1_lrg = std(GL_51_1_lrg.ev2o_ROC_AUC_1D_accum - GL_51_1_lrg.im_ROC_AUC_1D_accum);
        sGLev_51_1_mid = std(GL_51_1_mid.ev2o_ROC_AUC_1D_accum - GL_51_1_mid.im_ROC_AUC_1D_accum);
        sGLev_51_1_sml = std(GL_51_1_sml.ev2o_ROC_AUC_1D_accum - GL_51_1_sml.im_ROC_AUC_1D_accum);
        %
        sGLev_51_3_lrg = std(GL_51_3_lrg.ev2o_ROC_AUC_1D_accum - GL_51_3_lrg.im_ROC_AUC_1D_accum);
        sGLev_51_3_mid = std(GL_51_3_mid.ev2o_ROC_AUC_1D_accum - GL_51_3_mid.im_ROC_AUC_1D_accum);
        sGLev_51_3_sml = std(GL_51_3_sml.ev2o_ROC_AUC_1D_accum - GL_51_3_sml.im_ROC_AUC_1D_accum);
        %
        sGLev_51_5_lrg = std(GL_51_5_lrg.ev2o_ROC_AUC_1D_accum - GL_51_5_lrg.im_ROC_AUC_1D_accum);
        sGLev_51_5_mid = std(GL_51_5_mid.ev2o_ROC_AUC_1D_accum - GL_51_5_mid.im_ROC_AUC_1D_accum);
        sGLev_51_5_sml = std(GL_51_5_sml.ev2o_ROC_AUC_1D_accum - GL_51_5_sml.im_ROC_AUC_1D_accum);
        %
        sGLev_51_10_lrg = std(GL_51_10_lrg.ev2o_ROC_AUC_1D_accum - GL_51_10_lrg.im_ROC_AUC_1D_accum);
        sGLev_51_10_mid = std(GL_51_10_mid.ev2o_ROC_AUC_1D_accum - GL_51_10_mid.im_ROC_AUC_1D_accum);
        sGLev_51_10_sml = std(GL_51_10_sml.ev2o_ROC_AUC_1D_accum - GL_51_10_sml.im_ROC_AUC_1D_accum);


        % histograms for larger (101x101) image patches
        hGLev_101_1_lrg = hist(GL_101_1_lrg.ev2o_ROC_AUC_1D_accum - GL_101_1_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_101_1_mid = hist(GL_101_1_mid.ev2o_ROC_AUC_1D_accum - GL_101_1_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_101_1_sml = hist(GL_101_1_sml.ev2o_ROC_AUC_1D_accum - GL_101_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGLev_101_3_lrg = hist(GL_101_3_lrg.ev2o_ROC_AUC_1D_accum - GL_101_3_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_101_3_mid = hist(GL_101_3_mid.ev2o_ROC_AUC_1D_accum - GL_101_3_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_101_3_sml = hist(GL_101_3_sml.ev2o_ROC_AUC_1D_accum - GL_101_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGLev_101_5_lrg = hist(GL_101_5_lrg.ev2o_ROC_AUC_1D_accum - GL_101_5_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_101_5_mid = hist(GL_101_5_mid.ev2o_ROC_AUC_1D_accum - GL_101_5_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_101_5_sml = hist(GL_101_5_sml.ev2o_ROC_AUC_1D_accum - GL_101_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hGLev_101_10_lrg = hist(GL_101_10_lrg.ev2o_ROC_AUC_1D_accum - GL_101_10_lrg.im_ROC_AUC_1D_accum, bins);
        hGLev_101_10_mid = hist(GL_101_10_mid.ev2o_ROC_AUC_1D_accum - GL_101_10_mid.im_ROC_AUC_1D_accum, bins);
        hGLev_101_10_sml = hist(GL_101_10_sml.ev2o_ROC_AUC_1D_accum - GL_101_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (101x101) image patches
        mGLev_101_1_lrg = mean(GL_101_1_lrg.ev2o_ROC_AUC_1D_accum - GL_101_1_lrg.im_ROC_AUC_1D_accum);
        mGLev_101_1_mid = mean(GL_101_1_mid.ev2o_ROC_AUC_1D_accum - GL_101_1_mid.im_ROC_AUC_1D_accum);
        mGLev_101_1_sml = mean(GL_101_1_sml.ev2o_ROC_AUC_1D_accum - GL_101_1_sml.im_ROC_AUC_1D_accum);
        %
        mGLev_101_3_lrg = mean(GL_101_3_lrg.ev2o_ROC_AUC_1D_accum - GL_101_3_lrg.im_ROC_AUC_1D_accum);
        mGLev_101_3_mid = mean(GL_101_3_mid.ev2o_ROC_AUC_1D_accum - GL_101_3_mid.im_ROC_AUC_1D_accum);
        mGLev_101_3_sml = mean(GL_101_3_sml.ev2o_ROC_AUC_1D_accum - GL_101_3_sml.im_ROC_AUC_1D_accum);
        %
        mGLev_101_5_lrg = mean(GL_101_5_lrg.ev2o_ROC_AUC_1D_accum - GL_101_5_lrg.im_ROC_AUC_1D_accum);
        mGLev_101_5_mid = mean(GL_101_5_mid.ev2o_ROC_AUC_1D_accum - GL_101_5_mid.im_ROC_AUC_1D_accum);
        mGLev_101_5_sml = mean(GL_101_5_sml.ev2o_ROC_AUC_1D_accum - GL_101_5_sml.im_ROC_AUC_1D_accum);
        %
        mGLev_101_10_lrg = mean(GL_101_10_lrg.ev2o_ROC_AUC_1D_accum - GL_101_10_lrg.im_ROC_AUC_1D_accum);
        mGLev_101_10_mid = mean(GL_101_10_mid.ev2o_ROC_AUC_1D_accum - GL_101_10_mid.im_ROC_AUC_1D_accum);
        mGLev_101_10_sml = mean(GL_101_10_sml.ev2o_ROC_AUC_1D_accum - GL_101_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (101x101) image patches
        sGLev_101_1_lrg = std(GL_101_1_lrg.ev2o_ROC_AUC_1D_accum - GL_101_1_lrg.im_ROC_AUC_1D_accum);
        sGLev_101_1_mid = std(GL_101_1_mid.ev2o_ROC_AUC_1D_accum - GL_101_1_mid.im_ROC_AUC_1D_accum);
        sGLev_101_1_sml = std(GL_101_1_sml.ev2o_ROC_AUC_1D_accum - GL_101_1_sml.im_ROC_AUC_1D_accum);
        %
        sGLev_101_3_lrg = std(GL_101_3_lrg.ev2o_ROC_AUC_1D_accum - GL_101_3_lrg.im_ROC_AUC_1D_accum);
        sGLev_101_3_mid = std(GL_101_3_mid.ev2o_ROC_AUC_1D_accum - GL_101_3_mid.im_ROC_AUC_1D_accum);
        sGLev_101_3_sml = std(GL_101_3_sml.ev2o_ROC_AUC_1D_accum - GL_101_3_sml.im_ROC_AUC_1D_accum);
        %
        sGLev_101_5_lrg = std(GL_101_5_lrg.ev2o_ROC_AUC_1D_accum - GL_101_5_lrg.im_ROC_AUC_1D_accum);
        sGLev_101_5_mid = std(GL_101_5_mid.ev2o_ROC_AUC_1D_accum - GL_101_5_mid.im_ROC_AUC_1D_accum);
        sGLev_101_5_sml = std(GL_101_5_sml.ev2o_ROC_AUC_1D_accum - GL_101_5_sml.im_ROC_AUC_1D_accum);
        %
        sGLev_101_10_lrg = std(GL_101_10_lrg.ev2o_ROC_AUC_1D_accum - GL_101_10_lrg.im_ROC_AUC_1D_accum);
        sGLev_101_10_mid = std(GL_101_10_mid.ev2o_ROC_AUC_1D_accum - GL_101_10_mid.im_ROC_AUC_1D_accum);
        sGLev_101_10_sml = std(GL_101_10_sml.ev2o_ROC_AUC_1D_accum - GL_101_10_sml.im_ROC_AUC_1D_accum);



        %% Histogram, Mean & Std for Kur Upper Bound for Modularity SKH Methods

        % histograms for smaller (51x51) image patches
        hSK_51_1_lrg = hist(SK_51_1_lrg.kur_ROC_AUC_GS_accum - SK_51_1_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_51_1_mid = hist(SK_51_1_mid.kur_ROC_AUC_GS_accum - SK_51_1_mid.im_ROC_AUC_1D_accum, bins);
        hSK_51_1_sml = hist(SK_51_1_sml.kur_ROC_AUC_GS_accum - SK_51_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hSK_51_3_lrg = hist(SK_51_3_lrg.kur_ROC_AUC_GS_accum - SK_51_3_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_51_3_mid = hist(SK_51_3_mid.kur_ROC_AUC_GS_accum - SK_51_3_mid.im_ROC_AUC_1D_accum, bins);
        hSK_51_3_sml = hist(SK_51_3_sml.kur_ROC_AUC_GS_accum - SK_51_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hSK_51_5_lrg = hist(SK_51_5_lrg.kur_ROC_AUC_GS_accum - SK_51_5_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_51_5_mid = hist(SK_51_5_mid.kur_ROC_AUC_GS_accum - SK_51_5_mid.im_ROC_AUC_1D_accum, bins);
        hSK_51_5_sml = hist(SK_51_5_sml.kur_ROC_AUC_GS_accum - SK_51_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hSK_51_10_lrg = hist(SK_51_10_lrg.kur_ROC_AUC_GS_accum - SK_51_10_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_51_10_mid = hist(SK_51_10_mid.kur_ROC_AUC_GS_accum - SK_51_10_mid.im_ROC_AUC_1D_accum, bins);
        hSK_51_10_sml = hist(SK_51_10_sml.kur_ROC_AUC_GS_accum - SK_51_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (51x51) image patches
        mSK_51_1_lrg = mean(SK_51_1_lrg.kur_ROC_AUC_GS_accum - SK_51_1_lrg.im_ROC_AUC_1D_accum);
        mSK_51_1_mid = mean(SK_51_1_mid.kur_ROC_AUC_GS_accum - SK_51_1_mid.im_ROC_AUC_1D_accum);
        mSK_51_1_sml = mean(SK_51_1_sml.kur_ROC_AUC_GS_accum - SK_51_1_sml.im_ROC_AUC_1D_accum);
        %
        mSK_51_3_lrg = mean(SK_51_3_lrg.kur_ROC_AUC_GS_accum - SK_51_3_lrg.im_ROC_AUC_1D_accum);
        mSK_51_3_mid = mean(SK_51_3_mid.kur_ROC_AUC_GS_accum - SK_51_3_mid.im_ROC_AUC_1D_accum);
        mSK_51_3_sml = mean(SK_51_3_sml.kur_ROC_AUC_GS_accum - SK_51_3_sml.im_ROC_AUC_1D_accum);
        %
        mSK_51_5_lrg = mean(SK_51_5_lrg.kur_ROC_AUC_GS_accum - SK_51_5_lrg.im_ROC_AUC_1D_accum);
        mSK_51_5_mid = mean(SK_51_5_mid.kur_ROC_AUC_GS_accum - SK_51_5_mid.im_ROC_AUC_1D_accum);
        mSK_51_5_sml = mean(SK_51_5_sml.kur_ROC_AUC_GS_accum - SK_51_5_sml.im_ROC_AUC_1D_accum);
        %
        mSK_51_10_lrg = mean(SK_51_10_lrg.kur_ROC_AUC_GS_accum - SK_51_10_lrg.im_ROC_AUC_1D_accum);
        mSK_51_10_mid = mean(SK_51_10_mid.kur_ROC_AUC_GS_accum - SK_51_10_mid.im_ROC_AUC_1D_accum);
        mSK_51_10_sml = mean(SK_51_10_sml.kur_ROC_AUC_GS_accum - SK_51_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (51x51) image patches
        sSK_51_1_lrg = std(SK_51_1_lrg.kur_ROC_AUC_GS_accum - SK_51_1_lrg.im_ROC_AUC_1D_accum);
        sSK_51_1_mid = std(SK_51_1_mid.kur_ROC_AUC_GS_accum - SK_51_1_mid.im_ROC_AUC_1D_accum);
        sSK_51_1_sml = std(SK_51_1_sml.kur_ROC_AUC_GS_accum - SK_51_1_sml.im_ROC_AUC_1D_accum);
        %
        sSK_51_3_lrg = std(SK_51_3_lrg.kur_ROC_AUC_GS_accum - SK_51_3_lrg.im_ROC_AUC_1D_accum);
        sSK_51_3_mid = std(SK_51_3_mid.kur_ROC_AUC_GS_accum - SK_51_3_mid.im_ROC_AUC_1D_accum);
        sSK_51_3_sml = std(SK_51_3_sml.kur_ROC_AUC_GS_accum - SK_51_3_sml.im_ROC_AUC_1D_accum);
        %
        sSK_51_5_lrg = std(SK_51_5_lrg.kur_ROC_AUC_GS_accum - SK_51_5_lrg.im_ROC_AUC_1D_accum);
        sSK_51_5_mid = std(SK_51_5_mid.kur_ROC_AUC_GS_accum - SK_51_5_mid.im_ROC_AUC_1D_accum);
        sSK_51_5_sml = std(SK_51_5_sml.kur_ROC_AUC_GS_accum - SK_51_5_sml.im_ROC_AUC_1D_accum);
        %
        sSK_51_10_lrg = std(SK_51_10_lrg.kur_ROC_AUC_GS_accum - SK_51_10_lrg.im_ROC_AUC_1D_accum);
        sSK_51_10_mid = std(SK_51_10_mid.kur_ROC_AUC_GS_accum - SK_51_10_mid.im_ROC_AUC_1D_accum);
        sSK_51_10_sml = std(SK_51_10_sml.kur_ROC_AUC_GS_accum - SK_51_10_sml.im_ROC_AUC_1D_accum);


        % histograms for larger (101x101) image patches
        hSK_101_1_lrg = hist(SK_101_1_lrg.kur_ROC_AUC_GS_accum - SK_101_1_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_101_1_mid = hist(SK_101_1_mid.kur_ROC_AUC_GS_accum - SK_101_1_mid.im_ROC_AUC_1D_accum, bins);
        hSK_101_1_sml = hist(SK_101_1_sml.kur_ROC_AUC_GS_accum - SK_101_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hSK_101_3_lrg = hist(SK_101_3_lrg.kur_ROC_AUC_GS_accum - SK_101_3_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_101_3_mid = hist(SK_101_3_mid.kur_ROC_AUC_GS_accum - SK_101_3_mid.im_ROC_AUC_1D_accum, bins);
        hSK_101_3_sml = hist(SK_101_3_sml.kur_ROC_AUC_GS_accum - SK_101_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hSK_101_5_lrg = hist(SK_101_5_lrg.kur_ROC_AUC_GS_accum - SK_101_5_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_101_5_mid = hist(SK_101_5_mid.kur_ROC_AUC_GS_accum - SK_101_5_mid.im_ROC_AUC_1D_accum, bins);
        hSK_101_5_sml = hist(SK_101_5_sml.kur_ROC_AUC_GS_accum - SK_101_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hSK_101_10_lrg = hist(SK_101_10_lrg.kur_ROC_AUC_GS_accum - SK_101_10_lrg.im_ROC_AUC_1D_accum, bins);
        hSK_101_10_mid = hist(SK_101_10_mid.kur_ROC_AUC_GS_accum - SK_101_10_mid.im_ROC_AUC_1D_accum, bins);
        hSK_101_10_sml = hist(SK_101_10_sml.kur_ROC_AUC_GS_accum - SK_101_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (101x101) image patches
        mSK_101_1_lrg = mean(SK_101_1_lrg.kur_ROC_AUC_GS_accum - SK_101_1_lrg.im_ROC_AUC_1D_accum);
        mSK_101_1_mid = mean(SK_101_1_mid.kur_ROC_AUC_GS_accum - SK_101_1_mid.im_ROC_AUC_1D_accum);
        mSK_101_1_sml = mean(SK_101_1_sml.kur_ROC_AUC_GS_accum - SK_101_1_sml.im_ROC_AUC_1D_accum);
        %
        mSK_101_3_lrg = mean(SK_101_3_lrg.kur_ROC_AUC_GS_accum - SK_101_3_lrg.im_ROC_AUC_1D_accum);
        mSK_101_3_mid = mean(SK_101_3_mid.kur_ROC_AUC_GS_accum - SK_101_3_mid.im_ROC_AUC_1D_accum);
        mSK_101_3_sml = mean(SK_101_3_sml.kur_ROC_AUC_GS_accum - SK_101_3_sml.im_ROC_AUC_1D_accum);
        %
        mSK_101_5_lrg = mean(SK_101_5_lrg.kur_ROC_AUC_GS_accum - SK_101_5_lrg.im_ROC_AUC_1D_accum);
        mSK_101_5_mid = mean(SK_101_5_mid.kur_ROC_AUC_GS_accum - SK_101_5_mid.im_ROC_AUC_1D_accum);
        mSK_101_5_sml = mean(SK_101_5_sml.kur_ROC_AUC_GS_accum - SK_101_5_sml.im_ROC_AUC_1D_accum);
        %
        mSK_101_10_lrg = mean(SK_101_10_lrg.kur_ROC_AUC_GS_accum - SK_101_10_lrg.im_ROC_AUC_1D_accum);
        mSK_101_10_mid = mean(SK_101_10_mid.kur_ROC_AUC_GS_accum - SK_101_10_mid.im_ROC_AUC_1D_accum);
        mSK_101_10_sml = mean(SK_101_10_sml.kur_ROC_AUC_GS_accum - SK_101_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (101x101) image patches
        sSK_101_1_lrg = std(SK_101_1_lrg.kur_ROC_AUC_GS_accum - SK_101_1_lrg.im_ROC_AUC_1D_accum);
        sSK_101_1_mid = std(SK_101_1_mid.kur_ROC_AUC_GS_accum - SK_101_1_mid.im_ROC_AUC_1D_accum);
        sSK_101_1_sml = std(SK_101_1_sml.kur_ROC_AUC_GS_accum - SK_101_1_sml.im_ROC_AUC_1D_accum);
        %
        sSK_101_3_lrg = std(SK_101_3_lrg.kur_ROC_AUC_GS_accum - SK_101_3_lrg.im_ROC_AUC_1D_accum);
        sSK_101_3_mid = std(SK_101_3_mid.kur_ROC_AUC_GS_accum - SK_101_3_mid.im_ROC_AUC_1D_accum);
        sSK_101_3_sml = std(SK_101_3_sml.kur_ROC_AUC_GS_accum - SK_101_3_sml.im_ROC_AUC_1D_accum);
        %
        sSK_101_5_lrg = std(SK_101_5_lrg.kur_ROC_AUC_GS_accum - SK_101_5_lrg.im_ROC_AUC_1D_accum);
        sSK_101_5_mid = std(SK_101_5_mid.kur_ROC_AUC_GS_accum - SK_101_5_mid.im_ROC_AUC_1D_accum);
        sSK_101_5_sml = std(SK_101_5_sml.kur_ROC_AUC_GS_accum - SK_101_5_sml.im_ROC_AUC_1D_accum);
        %
        sSK_101_10_lrg = std(SK_101_10_lrg.kur_ROC_AUC_GS_accum - SK_101_10_lrg.im_ROC_AUC_1D_accum);
        sSK_101_10_mid = std(SK_101_10_mid.kur_ROC_AUC_GS_accum - SK_101_10_mid.im_ROC_AUC_1D_accum);
        sSK_101_10_sml = std(SK_101_10_sml.kur_ROC_AUC_GS_accum - SK_101_10_sml.im_ROC_AUC_1D_accum);





        %% Histogram, Mean & Std for Kur Upper Bound for Modularity N&G Methods

        % histograms for smaller (51x51) image patches
        hNG_51_1_lrg = hist(NG_51_1_lrg.kur_ROC_AUC_GS_accum - NG_51_1_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_51_1_mid = hist(NG_51_1_mid.kur_ROC_AUC_GS_accum - NG_51_1_mid.im_ROC_AUC_1D_accum, bins);
        hNG_51_1_sml = hist(NG_51_1_sml.kur_ROC_AUC_GS_accum - NG_51_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hNG_51_3_lrg = hist(NG_51_3_lrg.kur_ROC_AUC_GS_accum - NG_51_3_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_51_3_mid = hist(NG_51_3_mid.kur_ROC_AUC_GS_accum - NG_51_3_mid.im_ROC_AUC_1D_accum, bins);
        hNG_51_3_sml = hist(NG_51_3_sml.kur_ROC_AUC_GS_accum - NG_51_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hNG_51_5_lrg = hist(NG_51_5_lrg.kur_ROC_AUC_GS_accum - NG_51_5_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_51_5_mid = hist(NG_51_5_mid.kur_ROC_AUC_GS_accum - NG_51_5_mid.im_ROC_AUC_1D_accum, bins);
        hNG_51_5_sml = hist(NG_51_5_sml.kur_ROC_AUC_GS_accum - NG_51_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hNG_51_10_lrg = hist(NG_51_10_lrg.kur_ROC_AUC_GS_accum - NG_51_10_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_51_10_mid = hist(NG_51_10_mid.kur_ROC_AUC_GS_accum - NG_51_10_mid.im_ROC_AUC_1D_accum, bins);
        hNG_51_10_sml = hist(NG_51_10_sml.kur_ROC_AUC_GS_accum - NG_51_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (51x51) image patches
        mNG_51_1_lrg = mean(NG_51_1_lrg.kur_ROC_AUC_GS_accum - NG_51_1_lrg.im_ROC_AUC_1D_accum);
        mNG_51_1_mid = mean(NG_51_1_mid.kur_ROC_AUC_GS_accum - NG_51_1_mid.im_ROC_AUC_1D_accum);
        mNG_51_1_sml = mean(NG_51_1_sml.kur_ROC_AUC_GS_accum - NG_51_1_sml.im_ROC_AUC_1D_accum);
        %
        mNG_51_3_lrg = mean(NG_51_3_lrg.kur_ROC_AUC_GS_accum - NG_51_3_lrg.im_ROC_AUC_1D_accum);
        mNG_51_3_mid = mean(NG_51_3_mid.kur_ROC_AUC_GS_accum - NG_51_3_mid.im_ROC_AUC_1D_accum);
        mNG_51_3_sml = mean(NG_51_3_sml.kur_ROC_AUC_GS_accum - NG_51_3_sml.im_ROC_AUC_1D_accum);
        %
        mNG_51_5_lrg = mean(NG_51_5_lrg.kur_ROC_AUC_GS_accum - NG_51_5_lrg.im_ROC_AUC_1D_accum);
        mNG_51_5_mid = mean(NG_51_5_mid.kur_ROC_AUC_GS_accum - NG_51_5_mid.im_ROC_AUC_1D_accum);
        mNG_51_5_sml = mean(NG_51_5_sml.kur_ROC_AUC_GS_accum - NG_51_5_sml.im_ROC_AUC_1D_accum);
        %
        mNG_51_10_lrg = mean(NG_51_10_lrg.kur_ROC_AUC_GS_accum - NG_51_10_lrg.im_ROC_AUC_1D_accum);
        mNG_51_10_mid = mean(NG_51_10_mid.kur_ROC_AUC_GS_accum - NG_51_10_mid.im_ROC_AUC_1D_accum);
        mNG_51_10_sml = mean(NG_51_10_sml.kur_ROC_AUC_GS_accum - NG_51_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (51x51) image patches
        sNG_51_1_lrg = std(NG_51_1_lrg.kur_ROC_AUC_GS_accum - NG_51_1_lrg.im_ROC_AUC_1D_accum);
        sNG_51_1_mid = std(NG_51_1_mid.kur_ROC_AUC_GS_accum - NG_51_1_mid.im_ROC_AUC_1D_accum);
        sNG_51_1_sml = std(NG_51_1_sml.kur_ROC_AUC_GS_accum - NG_51_1_sml.im_ROC_AUC_1D_accum);
        %
        sNG_51_3_lrg = std(NG_51_3_lrg.kur_ROC_AUC_GS_accum - NG_51_3_lrg.im_ROC_AUC_1D_accum);
        sNG_51_3_mid = std(NG_51_3_mid.kur_ROC_AUC_GS_accum - NG_51_3_mid.im_ROC_AUC_1D_accum);
        sNG_51_3_sml = std(NG_51_3_sml.kur_ROC_AUC_GS_accum - NG_51_3_sml.im_ROC_AUC_1D_accum);
        %
        sNG_51_5_lrg = std(NG_51_5_lrg.kur_ROC_AUC_GS_accum - NG_51_5_lrg.im_ROC_AUC_1D_accum);
        sNG_51_5_mid = std(NG_51_5_mid.kur_ROC_AUC_GS_accum - NG_51_5_mid.im_ROC_AUC_1D_accum);
        sNG_51_5_sml = std(NG_51_5_sml.kur_ROC_AUC_GS_accum - NG_51_5_sml.im_ROC_AUC_1D_accum);
        %
        sNG_51_10_lrg = std(NG_51_10_lrg.kur_ROC_AUC_GS_accum - NG_51_10_lrg.im_ROC_AUC_1D_accum);
        sNG_51_10_mid = std(NG_51_10_mid.kur_ROC_AUC_GS_accum - NG_51_10_mid.im_ROC_AUC_1D_accum);
        sNG_51_10_sml = std(NG_51_10_sml.kur_ROC_AUC_GS_accum - NG_51_10_sml.im_ROC_AUC_1D_accum);


        % histograms for larger (101x101) image patches
        hNG_101_1_lrg = hist(NG_101_1_lrg.kur_ROC_AUC_GS_accum - NG_101_1_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_101_1_mid = hist(NG_101_1_mid.kur_ROC_AUC_GS_accum - NG_101_1_mid.im_ROC_AUC_1D_accum, bins);
        hNG_101_1_sml = hist(NG_101_1_sml.kur_ROC_AUC_GS_accum - NG_101_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hNG_101_3_lrg = hist(NG_101_3_lrg.kur_ROC_AUC_GS_accum - NG_101_3_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_101_3_mid = hist(NG_101_3_mid.kur_ROC_AUC_GS_accum - NG_101_3_mid.im_ROC_AUC_1D_accum, bins);
        hNG_101_3_sml = hist(NG_101_3_sml.kur_ROC_AUC_GS_accum - NG_101_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hNG_101_5_lrg = hist(NG_101_5_lrg.kur_ROC_AUC_GS_accum - NG_101_5_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_101_5_mid = hist(NG_101_5_mid.kur_ROC_AUC_GS_accum - NG_101_5_mid.im_ROC_AUC_1D_accum, bins);
        hNG_101_5_sml = hist(NG_101_5_sml.kur_ROC_AUC_GS_accum - NG_101_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hNG_101_10_lrg = hist(NG_101_10_lrg.kur_ROC_AUC_GS_accum - NG_101_10_lrg.im_ROC_AUC_1D_accum, bins);
        hNG_101_10_mid = hist(NG_101_10_mid.kur_ROC_AUC_GS_accum - NG_101_10_mid.im_ROC_AUC_1D_accum, bins);
        hNG_101_10_sml = hist(NG_101_10_sml.kur_ROC_AUC_GS_accum - NG_101_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (101x101) image patches
        mNG_101_1_lrg = mean(NG_101_1_lrg.kur_ROC_AUC_GS_accum - NG_101_1_lrg.im_ROC_AUC_1D_accum);
        mNG_101_1_mid = mean(NG_101_1_mid.kur_ROC_AUC_GS_accum - NG_101_1_mid.im_ROC_AUC_1D_accum);
        mNG_101_1_sml = mean(NG_101_1_sml.kur_ROC_AUC_GS_accum - NG_101_1_sml.im_ROC_AUC_1D_accum);
        %
        mNG_101_3_lrg = mean(NG_101_3_lrg.kur_ROC_AUC_GS_accum - NG_101_3_lrg.im_ROC_AUC_1D_accum);
        mNG_101_3_mid = mean(NG_101_3_mid.kur_ROC_AUC_GS_accum - NG_101_3_mid.im_ROC_AUC_1D_accum);
        mNG_101_3_sml = mean(NG_101_3_sml.kur_ROC_AUC_GS_accum - NG_101_3_sml.im_ROC_AUC_1D_accum);
        %
        mNG_101_5_lrg = mean(NG_101_5_lrg.kur_ROC_AUC_GS_accum - NG_101_5_lrg.im_ROC_AUC_1D_accum);
        mNG_101_5_mid = mean(NG_101_5_mid.kur_ROC_AUC_GS_accum - NG_101_5_mid.im_ROC_AUC_1D_accum);
        mNG_101_5_sml = mean(NG_101_5_sml.kur_ROC_AUC_GS_accum - NG_101_5_sml.im_ROC_AUC_1D_accum);
        %
        mNG_101_10_lrg = mean(NG_101_10_lrg.kur_ROC_AUC_GS_accum - NG_101_10_lrg.im_ROC_AUC_1D_accum);
        mNG_101_10_mid = mean(NG_101_10_mid.kur_ROC_AUC_GS_accum - NG_101_10_mid.im_ROC_AUC_1D_accum);
        mNG_101_10_sml = mean(NG_101_10_sml.kur_ROC_AUC_GS_accum - NG_101_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (101x101) image patches
        sNG_101_1_lrg = std(NG_101_1_lrg.kur_ROC_AUC_GS_accum - NG_101_1_lrg.im_ROC_AUC_1D_accum);
        sNG_101_1_mid = std(NG_101_1_mid.kur_ROC_AUC_GS_accum - NG_101_1_mid.im_ROC_AUC_1D_accum);
        sNG_101_1_sml = std(NG_101_1_sml.kur_ROC_AUC_GS_accum - NG_101_1_sml.im_ROC_AUC_1D_accum);
        %
        sNG_101_3_lrg = std(NG_101_3_lrg.kur_ROC_AUC_GS_accum - NG_101_3_lrg.im_ROC_AUC_1D_accum);
        sNG_101_3_mid = std(NG_101_3_mid.kur_ROC_AUC_GS_accum - NG_101_3_mid.im_ROC_AUC_1D_accum);
        sNG_101_3_sml = std(NG_101_3_sml.kur_ROC_AUC_GS_accum - NG_101_3_sml.im_ROC_AUC_1D_accum);
        %
        sNG_101_5_lrg = std(NG_101_5_lrg.kur_ROC_AUC_GS_accum - NG_101_5_lrg.im_ROC_AUC_1D_accum);
        sNG_101_5_mid = std(NG_101_5_mid.kur_ROC_AUC_GS_accum - NG_101_5_mid.im_ROC_AUC_1D_accum);
        sNG_101_5_sml = std(NG_101_5_sml.kur_ROC_AUC_GS_accum - NG_101_5_sml.im_ROC_AUC_1D_accum);
        %
        sNG_101_10_lrg = std(NG_101_10_lrg.kur_ROC_AUC_GS_accum - NG_101_10_lrg.im_ROC_AUC_1D_accum);
        sNG_101_10_mid = std(NG_101_10_mid.kur_ROC_AUC_GS_accum - NG_101_10_mid.im_ROC_AUC_1D_accum);
        sNG_101_10_sml = std(NG_101_10_sml.kur_ROC_AUC_GS_accum - NG_101_10_sml.im_ROC_AUC_1D_accum);




        %% Histogram, Mean & Std for Kur Upper Bound for Avg Assoc Methods

        % histograms for smaller (51x51) image patches
        hAA_51_1_lrg = hist(AA_51_1_lrg.kur_ROC_AUC_GS_accum - AA_51_1_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_51_1_mid = hist(AA_51_1_mid.kur_ROC_AUC_GS_accum - AA_51_1_mid.im_ROC_AUC_1D_accum, bins);
        hAA_51_1_sml = hist(AA_51_1_sml.kur_ROC_AUC_GS_accum - AA_51_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hAA_51_3_lrg = hist(AA_51_3_lrg.kur_ROC_AUC_GS_accum - AA_51_3_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_51_3_mid = hist(AA_51_3_mid.kur_ROC_AUC_GS_accum - AA_51_3_mid.im_ROC_AUC_1D_accum, bins);
        hAA_51_3_sml = hist(AA_51_3_sml.kur_ROC_AUC_GS_accum - AA_51_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hAA_51_5_lrg = hist(AA_51_5_lrg.kur_ROC_AUC_GS_accum - AA_51_5_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_51_5_mid = hist(AA_51_5_mid.kur_ROC_AUC_GS_accum - AA_51_5_mid.im_ROC_AUC_1D_accum, bins);
        hAA_51_5_sml = hist(AA_51_5_sml.kur_ROC_AUC_GS_accum - AA_51_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hAA_51_10_lrg = hist(AA_51_10_lrg.kur_ROC_AUC_GS_accum - AA_51_10_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_51_10_mid = hist(AA_51_10_mid.kur_ROC_AUC_GS_accum - AA_51_10_mid.im_ROC_AUC_1D_accum, bins);
        hAA_51_10_sml = hist(AA_51_10_sml.kur_ROC_AUC_GS_accum - AA_51_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (51x51) image patches
        mAA_51_1_lrg = mean(AA_51_1_lrg.kur_ROC_AUC_GS_accum - AA_51_1_lrg.im_ROC_AUC_1D_accum);
        mAA_51_1_mid = mean(AA_51_1_mid.kur_ROC_AUC_GS_accum - AA_51_1_mid.im_ROC_AUC_1D_accum);
        mAA_51_1_sml = mean(AA_51_1_sml.kur_ROC_AUC_GS_accum - AA_51_1_sml.im_ROC_AUC_1D_accum);
        %
        mAA_51_3_lrg = mean(AA_51_3_lrg.kur_ROC_AUC_GS_accum - AA_51_3_lrg.im_ROC_AUC_1D_accum);
        mAA_51_3_mid = mean(AA_51_3_mid.kur_ROC_AUC_GS_accum - AA_51_3_mid.im_ROC_AUC_1D_accum);
        mAA_51_3_sml = mean(AA_51_3_sml.kur_ROC_AUC_GS_accum - AA_51_3_sml.im_ROC_AUC_1D_accum);
        %
        mAA_51_5_lrg = mean(AA_51_5_lrg.kur_ROC_AUC_GS_accum - AA_51_5_lrg.im_ROC_AUC_1D_accum);
        mAA_51_5_mid = mean(AA_51_5_mid.kur_ROC_AUC_GS_accum - AA_51_5_mid.im_ROC_AUC_1D_accum);
        mAA_51_5_sml = mean(AA_51_5_sml.kur_ROC_AUC_GS_accum - AA_51_5_sml.im_ROC_AUC_1D_accum);
        %
        mAA_51_10_lrg = mean(AA_51_10_lrg.kur_ROC_AUC_GS_accum - AA_51_10_lrg.im_ROC_AUC_1D_accum);
        mAA_51_10_mid = mean(AA_51_10_mid.kur_ROC_AUC_GS_accum - AA_51_10_mid.im_ROC_AUC_1D_accum);
        mAA_51_10_sml = mean(AA_51_10_sml.kur_ROC_AUC_GS_accum - AA_51_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (51x51) image patches
        sAA_51_1_lrg = std(AA_51_1_lrg.kur_ROC_AUC_GS_accum - AA_51_1_lrg.im_ROC_AUC_1D_accum);
        sAA_51_1_mid = std(AA_51_1_mid.kur_ROC_AUC_GS_accum - AA_51_1_mid.im_ROC_AUC_1D_accum);
        sAA_51_1_sml = std(AA_51_1_sml.kur_ROC_AUC_GS_accum - AA_51_1_sml.im_ROC_AUC_1D_accum);
        %
        sAA_51_3_lrg = std(AA_51_3_lrg.kur_ROC_AUC_GS_accum - AA_51_3_lrg.im_ROC_AUC_1D_accum);
        sAA_51_3_mid = std(AA_51_3_mid.kur_ROC_AUC_GS_accum - AA_51_3_mid.im_ROC_AUC_1D_accum);
        sAA_51_3_sml = std(AA_51_3_sml.kur_ROC_AUC_GS_accum - AA_51_3_sml.im_ROC_AUC_1D_accum);
        %
        sAA_51_5_lrg = std(AA_51_5_lrg.kur_ROC_AUC_GS_accum - AA_51_5_lrg.im_ROC_AUC_1D_accum);
        sAA_51_5_mid = std(AA_51_5_mid.kur_ROC_AUC_GS_accum - AA_51_5_mid.im_ROC_AUC_1D_accum);
        sAA_51_5_sml = std(AA_51_5_sml.kur_ROC_AUC_GS_accum - AA_51_5_sml.im_ROC_AUC_1D_accum);
        %
        sAA_51_10_lrg = std(AA_51_10_lrg.kur_ROC_AUC_GS_accum - AA_51_10_lrg.im_ROC_AUC_1D_accum);
        sAA_51_10_mid = std(AA_51_10_mid.kur_ROC_AUC_GS_accum - AA_51_10_mid.im_ROC_AUC_1D_accum);
        sAA_51_10_sml = std(AA_51_10_sml.kur_ROC_AUC_GS_accum - AA_51_10_sml.im_ROC_AUC_1D_accum);


        % histograms for larger (101x101) image patches
        hAA_101_1_lrg = hist(AA_101_1_lrg.kur_ROC_AUC_GS_accum - AA_101_1_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_101_1_mid = hist(AA_101_1_mid.kur_ROC_AUC_GS_accum - AA_101_1_mid.im_ROC_AUC_1D_accum, bins);
        hAA_101_1_sml = hist(AA_101_1_sml.kur_ROC_AUC_GS_accum - AA_101_1_sml.im_ROC_AUC_1D_accum, bins);
        %
        hAA_101_3_lrg = hist(AA_101_3_lrg.kur_ROC_AUC_GS_accum - AA_101_3_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_101_3_mid = hist(AA_101_3_mid.kur_ROC_AUC_GS_accum - AA_101_3_mid.im_ROC_AUC_1D_accum, bins);
        hAA_101_3_sml = hist(AA_101_3_sml.kur_ROC_AUC_GS_accum - AA_101_3_sml.im_ROC_AUC_1D_accum, bins);
        %
        hAA_101_5_lrg = hist(AA_101_5_lrg.kur_ROC_AUC_GS_accum - AA_101_5_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_101_5_mid = hist(AA_101_5_mid.kur_ROC_AUC_GS_accum - AA_101_5_mid.im_ROC_AUC_1D_accum, bins);
        hAA_101_5_sml = hist(AA_101_5_sml.kur_ROC_AUC_GS_accum - AA_101_5_sml.im_ROC_AUC_1D_accum, bins);
        %
        hAA_101_10_lrg = hist(AA_101_10_lrg.kur_ROC_AUC_GS_accum - AA_101_10_lrg.im_ROC_AUC_1D_accum, bins);
        hAA_101_10_mid = hist(AA_101_10_mid.kur_ROC_AUC_GS_accum - AA_101_10_mid.im_ROC_AUC_1D_accum, bins);
        hAA_101_10_sml = hist(AA_101_10_sml.kur_ROC_AUC_GS_accum - AA_101_10_sml.im_ROC_AUC_1D_accum, bins);


        % mean of distributions of smaller (101x101) image patches
        mAA_101_1_lrg = mean(AA_101_1_lrg.kur_ROC_AUC_GS_accum - AA_101_1_lrg.im_ROC_AUC_1D_accum);
        mAA_101_1_mid = mean(AA_101_1_mid.kur_ROC_AUC_GS_accum - AA_101_1_mid.im_ROC_AUC_1D_accum);
        mAA_101_1_sml = mean(AA_101_1_sml.kur_ROC_AUC_GS_accum - AA_101_1_sml.im_ROC_AUC_1D_accum);
        %
        mAA_101_3_lrg = mean(AA_101_3_lrg.kur_ROC_AUC_GS_accum - AA_101_3_lrg.im_ROC_AUC_1D_accum);
        mAA_101_3_mid = mean(AA_101_3_mid.kur_ROC_AUC_GS_accum - AA_101_3_mid.im_ROC_AUC_1D_accum);
        mAA_101_3_sml = mean(AA_101_3_sml.kur_ROC_AUC_GS_accum - AA_101_3_sml.im_ROC_AUC_1D_accum);
        %
        mAA_101_5_lrg = mean(AA_101_5_lrg.kur_ROC_AUC_GS_accum - AA_101_5_lrg.im_ROC_AUC_1D_accum);
        mAA_101_5_mid = mean(AA_101_5_mid.kur_ROC_AUC_GS_accum - AA_101_5_mid.im_ROC_AUC_1D_accum);
        mAA_101_5_sml = mean(AA_101_5_sml.kur_ROC_AUC_GS_accum - AA_101_5_sml.im_ROC_AUC_1D_accum);
        %
        mAA_101_10_lrg = mean(AA_101_10_lrg.kur_ROC_AUC_GS_accum - AA_101_10_lrg.im_ROC_AUC_1D_accum);
        mAA_101_10_mid = mean(AA_101_10_mid.kur_ROC_AUC_GS_accum - AA_101_10_mid.im_ROC_AUC_1D_accum);
        mAA_101_10_sml = mean(AA_101_10_sml.kur_ROC_AUC_GS_accum - AA_101_10_sml.im_ROC_AUC_1D_accum);

        % std of distributions of smaller (101x101) image patches
        sAA_101_1_lrg = std(AA_101_1_lrg.kur_ROC_AUC_GS_accum - AA_101_1_lrg.im_ROC_AUC_1D_accum);
        sAA_101_1_mid = std(AA_101_1_mid.kur_ROC_AUC_GS_accum - AA_101_1_mid.im_ROC_AUC_1D_accum);
        sAA_101_1_sml = std(AA_101_1_sml.kur_ROC_AUC_GS_accum - AA_101_1_sml.im_ROC_AUC_1D_accum);
        %
        sAA_101_3_lrg = std(AA_101_3_lrg.kur_ROC_AUC_GS_accum - AA_101_3_lrg.im_ROC_AUC_1D_accum);
        sAA_101_3_mid = std(AA_101_3_mid.kur_ROC_AUC_GS_accum - AA_101_3_mid.im_ROC_AUC_1D_accum);
        sAA_101_3_sml = std(AA_101_3_sml.kur_ROC_AUC_GS_accum - AA_101_3_sml.im_ROC_AUC_1D_accum);
        %
        sAA_101_5_lrg = std(AA_101_5_lrg.kur_ROC_AUC_GS_accum - AA_101_5_lrg.im_ROC_AUC_1D_accum);
        sAA_101_5_mid = std(AA_101_5_mid.kur_ROC_AUC_GS_accum - AA_101_5_mid.im_ROC_AUC_1D_accum);
        sAA_101_5_sml = std(AA_101_5_sml.kur_ROC_AUC_GS_accum - AA_101_5_sml.im_ROC_AUC_1D_accum);
        %
        sAA_101_10_lrg = std(AA_101_10_lrg.kur_ROC_AUC_GS_accum - AA_101_10_lrg.im_ROC_AUC_1D_accum);
        sAA_101_10_mid = std(AA_101_10_mid.kur_ROC_AUC_GS_accum - AA_101_10_mid.im_ROC_AUC_1D_accum);
        sAA_101_10_sml = std(AA_101_10_sml.kur_ROC_AUC_GS_accum - AA_101_10_sml.im_ROC_AUC_1D_accum);



        %% Plot mean & std of Isotropic Diffusion Kuramoto performance for rM & KS values to determine best performer.
        if(0)
            figure, hold on,
            errorbar(0.0, mIso_101_1_sml, sIso_101_1_sml, 'b.')
            errorbar(0.2, mIso_101_1_mid, sIso_101_1_mid, 'g.')
            errorbar(0.4, mIso_101_1_lrg, sIso_101_1_lrg, 'r.')
            errorbar(0.1, mIso_51_1_sml, sIso_51_1_sml, 'bo')
            errorbar(0.3, mIso_51_1_mid, sIso_51_1_mid, 'go')
            errorbar(0.5, mIso_51_1_lrg, sIso_51_1_lrg, 'ro')
            %
            %
            errorbar(1.0, mIso_101_3_sml, sIso_101_3_sml, 'b.')
            errorbar(1.2, mIso_101_3_mid, sIso_101_3_mid, 'g.')
            errorbar(1.4, mIso_101_3_lrg, sIso_101_3_lrg, 'r.')
            errorbar(1.1, mIso_51_3_sml, sIso_51_3_sml, 'bo')
            errorbar(1.3, mIso_51_3_mid, sIso_51_3_mid, 'go')
            errorbar(1.5, mIso_51_3_lrg, sIso_51_3_lrg, 'ro')
            %
            errorbar(2.0, mIso_101_5_sml, sIso_101_5_sml, 'b.')
            errorbar(2.2, mIso_101_5_mid, sIso_101_5_mid, 'g.')
            errorbar(2.4, mIso_101_5_lrg, sIso_101_5_lrg, 'r.')
            errorbar(2.1, mIso_51_5_sml, sIso_51_5_sml, 'bo')
            errorbar(2.3, mIso_51_5_mid, sIso_51_5_mid, 'go')
            errorbar(2.5, mIso_51_5_lrg, sIso_51_5_lrg, 'ro')
            %
            errorbar(3.0, mIso_101_10_sml, sIso_101_10_sml, 'b.')
            errorbar(3.2, mIso_101_10_mid, sIso_101_10_mid, 'g.')
            errorbar(3.4, mIso_101_10_lrg, sIso_101_10_lrg, 'r.')
            errorbar(3.1, mIso_51_10_sml, sIso_51_10_sml, 'bo')
            errorbar(3.3, mIso_51_10_mid, sIso_51_10_mid, 'go')
            errorbar(3.5, mIso_51_10_lrg, sIso_51_10_lrg, 'ro')
            %
            legend({'ks sml','ks mid','ks lrg'})
            set(gca,'XTick',[0,1,2,3],'XTickLabel',[1,3,5,10],'FontSize',16,'FontWeight','Bold')
            xlabel('rM value','FontSize',18,'FontWeight','Bold')
            ylabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            title('(1). What are best params for Isotropic Diffusion Kuramoto? (large img patch)','FontSize',18,'FontWeight','Bold')
            grid on
            plot([0 3.5],[0 0],'k--')
        end





        %% Plot Isotropic Diffusion Kuramoto Histogram results for different img patch size and different KS
        if(0)
            figure, 
            %
            subplot(211),hold on
            plot(bins, normlze(hIso_101_3_lrg),'r')
            plot(bins, normlze(hIso_101_3_mid),'g')
            plot(bins, normlze(hIso_101_3_sml),'b')
            %
            plot(bins, normlze(hIso_51_3_lrg) ,'r--')
            plot(bins, normlze(hIso_51_3_mid) ,'g--')
            plot(bins, normlze(hIso_51_3_sml) ,'b--')
            %
            herrorbar(mIso_101_3_lrg, 0.32, sIso_101_3_lrg, 'r.')
            herrorbar(mIso_101_3_mid, 0.325, sIso_101_3_mid, 'g.')
            herrorbar(mIso_101_3_sml, 0.33, sIso_101_3_sml, 'b.')
            %
            plot([0 0], [0 0.3], 'k--','LineWidth',1.2)
            title('(rM=3)','FontSize',16,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            %
            subplot(212), hold on
            plot(bins, normlze(hIso_101_10_lrg),'r')
            plot(bins, normlze(hIso_101_10_mid),'g')
            plot(bins, normlze(hIso_101_10_sml),'b')
            %
            plot(bins, normlze(hIso_51_10_lrg) ,'r--')
            plot(bins, normlze(hIso_51_10_mid) ,'g--')
            plot(bins, normlze(hIso_51_10_sml) ,'b--')
            legend({'lrg ks lrg img','mid ks lrg img','sml ks lrg img','lrg ks sml img','mid ks sml img','sml ks sml img'})
            %
            herrorbar(mIso_101_10_lrg, 0.32, sIso_101_10_lrg, 'r.')
            herrorbar(mIso_101_10_mid, 0.325, sIso_101_10_mid, 'g.')
            herrorbar(mIso_101_10_sml, 0.33, sIso_101_10_sml, 'b.')
            %
            plot([0 0], [0 0.3], 'k--','LineWidth',1.2)
            title('(rM=10)','FontSize',16,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            xlabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold') 
            ylabel('counts (segment pairs)','FontSize',18,'FontWeight','Bold')
            %
            annotation('textbox', [0 0.9 1 0.1],'String',['(2). Isotropic Diffusion Kuramoto : How does ks scale & Img patch size effect Segmentation?'], ...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
        end





        %% Plot mean & std of Graph Laplacian Kuramoto performance for rM & KS values to determine best performer.
        if(0)
            figure, hold on,
            errorbar(0.0, mGL_101_1_sml, sGL_101_1_sml, 'b.')
            errorbar(0.2, mGL_101_1_mid, sGL_101_1_mid, 'g.')
            errorbar(0.4, mGL_101_1_lrg, sGL_101_1_lrg, 'r.')
            errorbar(0.1, mGL_51_1_sml, sGL_51_1_sml, 'bo')
            errorbar(0.3, mGL_51_1_mid, sGL_51_1_mid, 'go')
            errorbar(0.5, mGL_51_1_lrg, sGL_51_1_lrg, 'ro')
            %
            %
            errorbar(1.0, mGL_101_3_sml, sGL_101_3_sml, 'b.')
            errorbar(1.2, mGL_101_3_mid, sGL_101_3_mid, 'g.')
            errorbar(1.4, mGL_101_3_lrg, sGL_101_3_lrg, 'r.')
            errorbar(1.1, mGL_51_3_sml, sGL_51_3_sml, 'bo')
            errorbar(1.3, mGL_51_3_mid, sGL_51_3_mid, 'go')
            errorbar(1.5, mGL_51_3_lrg, sGL_51_3_lrg, 'ro')
            %
            errorbar(2.0, mGL_101_5_sml, sGL_101_5_sml, 'b.')
            errorbar(2.2, mGL_101_5_mid, sGL_101_5_mid, 'g.')
            errorbar(2.4, mGL_101_5_lrg, sGL_101_5_lrg, 'r.')
            errorbar(2.1, mGL_51_5_sml, sGL_51_5_sml, 'bo')
            errorbar(2.3, mGL_51_5_mid, sGL_51_5_mid, 'go')
            errorbar(2.5, mGL_51_5_lrg, sGL_51_5_lrg, 'ro')
            %
            errorbar(3.0, mGL_101_10_sml, sGL_101_10_sml, 'b.')
            errorbar(3.2, mGL_101_10_mid, sGL_101_10_mid, 'g.')
            errorbar(3.4, mGL_101_10_lrg, sGL_101_10_lrg, 'r.')
            errorbar(3.1, mGL_51_10_sml, sGL_51_10_sml, 'bo')
            errorbar(3.3, mGL_51_10_mid, sGL_51_10_mid, 'go')
            errorbar(3.5, mGL_51_10_lrg, sGL_51_10_lrg, 'ro')
            %
            legend({'ks sml','ks mid','ks lrg'})
            set(gca,'XTick',[0,1,2,3],'XTickLabel',[1,3,5,10],'FontSize',16,'FontWeight','Bold')
            xlabel('rM value','FontSize',18,'FontWeight','Bold')
            ylabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            title('(1). What are best params for Graph Laplacain Kuramoto?','FontSize',18,'FontWeight','Bold')
            grid on
            plot([0 3.5],[0 0],'k--')
        end


        %% Plot mean & std of Graph Laplacian EigVec2 performance for rM & KS values to determine best performer.  
        if(0)
            figure, hold on,
            errorbar(0.0, mGLev_101_1_sml, sGLev_101_1_sml, 'b.')
            errorbar(0.2, mGLev_101_1_mid, sGLev_101_1_mid, 'g.')
            errorbar(0.4, mGLev_101_1_lrg, sGLev_101_1_lrg, 'r.')
            errorbar(0.1, mGLev_51_1_sml, sGLev_51_1_sml, 'bo')
            errorbar(0.3, mGLev_51_1_mid, sGLev_51_1_mid, 'go')
            errorbar(0.5, mGLev_51_1_lrg, sGLev_51_1_lrg, 'ro')
            %
            %
            errorbar(1.0, mGLev_101_3_sml, sGLev_101_3_sml, 'b.')
            errorbar(1.2, mGLev_101_3_mid, sGLev_101_3_mid, 'g.')
            errorbar(1.4, mGLev_101_3_lrg, sGLev_101_3_lrg, 'r.')
            errorbar(1.1, mGLev_51_3_sml, sGLev_51_3_sml, 'bo')
            errorbar(1.3, mGLev_51_3_mid, sGLev_51_3_mid, 'go')
            errorbar(1.5, mGLev_51_3_lrg, sGLev_51_3_lrg, 'ro')
            %
            errorbar(2.0, mGLev_101_5_sml, sGLev_101_5_sml, 'b.')
            errorbar(2.2, mGLev_101_5_mid, sGLev_101_5_mid, 'g.')
            errorbar(2.4, mGLev_101_5_lrg, sGLev_101_5_lrg, 'r.')
            errorbar(2.1, mGLev_51_5_sml, sGLev_51_5_sml, 'bo')
            errorbar(2.3, mGLev_51_5_mid, sGLev_51_5_mid, 'go')
            errorbar(2.5, mGLev_51_5_lrg, sGLev_51_5_lrg, 'ro')
            %
            errorbar(3.0, mGLev_101_10_sml, sGLev_101_10_sml, 'b.')
            errorbar(3.2, mGLev_101_10_mid, sGLev_101_10_mid, 'g.')
            errorbar(3.4, mGLev_101_10_lrg, sGLev_101_10_lrg, 'r.')
            errorbar(3.1, mGLev_51_10_sml, sGLev_51_10_sml, 'bo')
            errorbar(3.3, mGLev_51_10_mid, sGLev_51_10_mid, 'go')
            errorbar(3.5, mGLev_51_10_lrg, sGLev_51_10_lrg, 'ro')
            %
            legend({'ks sml','ks mid','ks lrg'})
            set(gca,'XTick',[0,1,2,3],'XTickLabel',[1,3,5,10],'FontSize',16,'FontWeight','Bold')
            xlabel('rM value','FontSize',18,'FontWeight','Bold')
            ylabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            title('(1). What are best params for Graph Laplacain EigVec2?','FontSize',18,'FontWeight','Bold')
            grid on
            plot([0 3.5],[0 0],'k--')
        end


        %% Plot mean & std of Modularity SK Kuramoto performance for rM & KS values to determine best performer.  
        if(0)
            figure, hold on,
            errorbar(0.0, mSK_101_1_sml, sSK_101_1_sml, 'b.')
            errorbar(0.2, mSK_101_1_mid, sSK_101_1_mid, 'g.')
            errorbar(0.4, mSK_101_1_lrg, sSK_101_1_lrg, 'r.')
            errorbar(0.1, mSK_51_1_sml, sSK_51_1_sml, 'bo')
            errorbar(0.3, mSK_51_1_mid, sSK_51_1_mid, 'go')
            errorbar(0.5, mSK_51_1_lrg, sSK_51_1_lrg, 'ro')
            %
            %
            errorbar(1.0, mSK_101_3_sml, sSK_101_3_sml, 'b.')
            errorbar(1.2, mSK_101_3_mid, sSK_101_3_mid, 'g.')
            errorbar(1.4, mSK_101_3_lrg, sSK_101_3_lrg, 'r.')
            errorbar(1.1, mSK_51_3_sml, sSK_51_3_sml, 'bo')
            errorbar(1.3, mSK_51_3_mid, sSK_51_3_mid, 'go')
            errorbar(1.5, mSK_51_3_lrg, sSK_51_3_lrg, 'ro')
            %
            errorbar(2.0, mSK_101_5_sml, sSK_101_5_sml, 'b.')
            errorbar(2.2, mSK_101_5_mid, sSK_101_5_mid, 'g.')
            errorbar(2.4, mSK_101_5_lrg, sSK_101_5_lrg, 'r.')
            errorbar(2.1, mSK_51_5_sml, sSK_51_5_sml, 'bo')
            errorbar(2.3, mSK_51_5_mid, sSK_51_5_mid, 'go')
            errorbar(2.5, mSK_51_5_lrg, sSK_51_5_lrg, 'ro')
            %
            errorbar(3.0, mSK_101_10_sml, sSK_101_10_sml, 'b.')
            errorbar(3.2, mSK_101_10_mid, sSK_101_10_mid, 'g.')
            errorbar(3.4, mSK_101_10_lrg, sSK_101_10_lrg, 'r.')
            errorbar(3.1, mSK_51_10_sml, sSK_51_10_sml, 'bo')
            errorbar(3.3, mSK_51_10_mid, sSK_51_10_mid, 'go')
            errorbar(3.5, mSK_51_10_lrg, sSK_51_10_lrg, 'ro')
            %
            legend({'ks sml','ks mid','ks lrg'})
            set(gca,'XTick',[0,1,2,3],'XTickLabel',[1,3,5,10],'FontSize',16,'FontWeight','Bold')
            xlabel('rM value','FontSize',18,'FontWeight','Bold')
            ylabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            title('(1). What are best params for Modularity SKH Kuramoto?','FontSize',18,'FontWeight','Bold')
            grid on
            plot([0 3.5],[0 0],'k--')
        end


        %% Plot mean & std of Modularity NG Kuramoto performance for rM & KS values to determine best performer.  
        if(0)
            figure, hold on,
            errorbar(0.0, mNG_101_1_sml, sNG_101_1_sml, 'b.')
            errorbar(0.2, mNG_101_1_mid, sNG_101_1_mid, 'g.')
            errorbar(0.4, mNG_101_1_lrg, sNG_101_1_lrg, 'r.')
            errorbar(0.1, mNG_51_1_sml, sNG_51_1_sml, 'bo')
            errorbar(0.3, mNG_51_1_mid, sNG_51_1_mid, 'go')
            errorbar(0.5, mNG_51_1_lrg, sNG_51_1_lrg, 'ro')
            %
            %
            errorbar(1.0, mNG_101_3_sml, sNG_101_3_sml, 'b.')
            errorbar(1.2, mNG_101_3_mid, sNG_101_3_mid, 'g.')
            errorbar(1.4, mNG_101_3_lrg, sNG_101_3_lrg, 'r.')
            errorbar(1.1, mNG_51_3_sml, sNG_51_3_sml, 'bo')
            errorbar(1.3, mNG_51_3_mid, sNG_51_3_mid, 'go')
            errorbar(1.5, mNG_51_3_lrg, sNG_51_3_lrg, 'ro')
            %
            errorbar(2.0, mNG_101_5_sml, sNG_101_5_sml, 'b.')
            errorbar(2.2, mNG_101_5_mid, sNG_101_5_mid, 'g.')
            errorbar(2.4, mNG_101_5_lrg, sNG_101_5_lrg, 'r.')
            errorbar(2.1, mNG_51_5_sml, sNG_51_5_sml, 'bo')
            errorbar(2.3, mNG_51_5_mid, sNG_51_5_mid, 'go')
            errorbar(2.5, mNG_51_5_lrg, sNG_51_5_lrg, 'ro')
            %
            errorbar(3.0, mNG_101_10_sml, sNG_101_10_sml, 'b.')
            errorbar(3.2, mNG_101_10_mid, sNG_101_10_mid, 'g.')
            errorbar(3.4, mNG_101_10_lrg, sNG_101_10_lrg, 'r.')
            errorbar(3.1, mNG_51_10_sml, sNG_51_10_sml, 'bo')
            errorbar(3.3, mNG_51_10_mid, sNG_51_10_mid, 'go')
            errorbar(3.5, mNG_51_10_lrg, sNG_51_10_lrg, 'ro')
            %
            legend({'ks sml','ks mid','ks lrg'})
            set(gca,'XTick',[0,1,2,3],'XTickLabel',[1,3,5,10],'FontSize',16,'FontWeight','Bold')
            xlabel('rM value','FontSize',18,'FontWeight','Bold')
            ylabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            title('(1). What are best params for Modularity N&G Kuramoto?','FontSize',18,'FontWeight','Bold')
            grid on
            plot([0 3.5],[0 0],'k--')
        end




        %% Plot mean & std of Avg Association Kuramoto performance for rM & KS values to determine best performer.  
        if(0)
            figure, hold on,
            errorbar(0.0, mAA_101_1_sml, sAA_101_1_sml, 'b.')
            errorbar(0.2, mAA_101_1_mid, sAA_101_1_mid, 'g.')
            errorbar(0.4, mAA_101_1_lrg, sAA_101_1_lrg, 'r.')
            errorbar(0.1, mAA_51_1_sml, sAA_51_1_sml, 'bo')
            errorbar(0.3, mAA_51_1_mid, sAA_51_1_mid, 'go')
            errorbar(0.5, mAA_51_1_lrg, sAA_51_1_lrg, 'ro')
            %
            %
            errorbar(1.0, mAA_101_3_sml, sAA_101_3_sml, 'b.')
            errorbar(1.2, mAA_101_3_mid, sAA_101_3_mid, 'g.')
            errorbar(1.4, mAA_101_3_lrg, sAA_101_3_lrg, 'r.')
            errorbar(1.1, mAA_51_3_sml, sAA_51_3_sml, 'bo')
            errorbar(1.3, mAA_51_3_mid, sAA_51_3_mid, 'go')
            errorbar(1.5, mAA_51_3_lrg, sAA_51_3_lrg, 'ro')
            %
            errorbar(2.0, mAA_101_5_sml, sAA_101_5_sml, 'b.')
            errorbar(2.2, mAA_101_5_mid, sAA_101_5_mid, 'g.')
            errorbar(2.4, mAA_101_5_lrg, sAA_101_5_lrg, 'r.')
            errorbar(2.1, mAA_51_5_sml, sAA_51_5_sml, 'bo')
            errorbar(2.3, mAA_51_5_mid, sAA_51_5_mid, 'go')
            errorbar(2.5, mAA_51_5_lrg, sAA_51_5_lrg, 'ro')
            %
            errorbar(3.0, mAA_101_10_sml, sAA_101_10_sml, 'b.')
            errorbar(3.2, mAA_101_10_mid, sAA_101_10_mid, 'g.')
            errorbar(3.4, mAA_101_10_lrg, sAA_101_10_lrg, 'r.')
            errorbar(3.1, mAA_51_10_sml, sAA_51_10_sml, 'bo')
            errorbar(3.3, mAA_51_10_mid, sAA_51_10_mid, 'go')
            errorbar(3.5, mAA_51_10_lrg, sAA_51_10_lrg, 'ro')
            %
            legend({'ks sml','ks mid','ks lrg'})
            set(gca,'XTick',[0,1,2,3],'XTickLabel',[1,3,5,10],'FontSize',16,'FontWeight','Bold')
            xlabel('rM value','FontSize',18,'FontWeight','Bold')
            ylabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            title('(1). What are best params for Avg Assoc Kuramoto?','FontSize',18,'FontWeight','Bold')
            grid on
            plot([0 3.5],[0 0],'k--')
        end



        %% Plot Best results Histograms for all Kuramoto methods
        figure, 
        subplot(211),hold on, 
        plot(bins, normlze(hIso_101_3_mid),'g','LineWidth',2)       % IsoDiff Kuramoto
        plot(bins, normlze(hSK_101_5_lrg),'r','LineWidth',2)         % SK Kuramoto
        plot(bins, normlze(hNG_101_10_lrg),'c','LineWidth',2)       % NG Kuramoto
        plot(bins, normlze(hGL_101_1_lrg),'m','LineWidth',2)         % GL Kuramoto
        %plot(bins, normlze(hGLev_101_3_lrg),'m--','LineWidth',2)  % GL 2nd eigenvector
        plot(bins, normlze(hAA_101_1_lrg),'y','LineWidth',2)         % AA Kuramoto
        %
        plot([0 0], [0 0.4], 'k--','LineWidth',1.1)
        xlabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
        ylabel('counts (segment pairs)','FontSize',18,'FontWeight','Bold')
        title('Best Performance','FontSize',16,'FontWeight','Bold')
        set(gca,'FontSize',16,'FontWeight','Bold')
        legend({'IsoDiff Kur','SK Kur','NG Kur','GL Kur','AA Kur'}) % 'GL Ev2'
        xlim([-0.6 0.6])
        %
        herrorbar(mIso_101_3_mid, 0.47, sIso_101_3_mid, 'g.') 
        herrorbar(mSK_101_5_lrg, 0.45, sSK_101_5_lrg, 'r.')
        herrorbar(mNG_101_10_lrg, 0.44, sNG_101_10_lrg, 'c.')
        herrorbar(mGL_101_3_lrg, 0.46, sGL_101_3_lrg, 'm.')
        herrorbar(mAA_101_1_lrg, 0.43, sAA_101_1_lrg, 'y.')
        text(0.1,0.48,'Mean & STD','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
        %
        %
        % Plot difference in mass distributions for Isotropic Diffusion and other Methods.
        subplot(212), 
        hold on,
        plot(bins, (hSK_101_5_lrg - hIso_101_3_mid) ./ sum(hIso_101_3_mid),'r','LineWidth',2)
        plot(bins, (hNG_101_10_lrg - hIso_101_3_mid) ./ sum(hIso_101_3_mid),'c','LineWidth',2)
        plot(bins, (hGL_101_1_lrg - hIso_101_3_mid) ./ sum(hIso_101_3_mid),'m','LineWidth',2)
        plot(bins, (hAA_101_1_lrg - hIso_101_3_mid) ./ sum(hIso_101_3_mid),'y','LineWidth',2)


        plot([-1 +1],[0 0],'k--')
        plot([0 0],[-0.05 +0.2], 'k--')
        title('Migration of mass from IsoDiff?','FontSize',18,'FontWeight','Bold')
        %legend({'SK Kur','NG Kur','GL Kur','AA Kur'})
        ylabel({'Mass Lost (-) and Gained (+)','from Isotropic Diffusion Distribution'},'FontSize',16,'FontWeight','Bold')
        xlabel('Performance Worse than (-) and Better than (+) just using Image Pixels','FontSize',16,'FontWeight','Bold')
        set(gca,'FontSize',16,'FontWeight','Bold')
        xlim([-0.5 +0.5])



        disp('% of points below diagonal line (worse than just image pixels)')
        disp('Iso:')
        sum(hIso_101_3_mid(1:10)) ./ sum(hIso_101_3_mid)
        disp('GL:')
        sum(hGL_101_1_lrg(1:10)) ./ sum(hGL_101_1_lrg)
        disp('SK:')
        sum(hSK_101_5_lrg(1:10)) ./ sum(hSK_101_5_lrg)
        disp('NG:')
        sum(hNG_101_10_lrg(1:10)) ./ sum(hNG_101_10_lrg)
        disp('AA:')
        sum(hAA_101_1_lrg(1:10)) ./ sum(hAA_101_1_lrg)







        %% Characterize how size of cluster pair relates to AUCm-AUCi metric

        % Iso Diff Distribution split by large and small cluster pairs
        if(0)
            IsoAUCdiff = Iso_101_3_mid.kur_ROC_AUC_GS_accum - Iso_101_3_mid.im_ROC_AUC_1D_accum;

            th_numpix = 1000;
            ind = (Iso_101_3_mid.ClusterSizes(1,:) < th_numpix & Iso_101_3_mid.ClusterSizes(2,:) < th_numpix);
            n = numel(find(ind>0));
            h1 = hist(IsoAUCdiff(ind),bins);
            h2 = hist(IsoAUCdiff(~ind),bins);

            figure, hold on,
            bar(bins,h2,0.50,'b')
            bar(bins,h1,0.25,'r')
            xlim([-0.55,0.55])
            legend({['Both Clusters > ',num2str(th_numpix),' pix'],['Either Cluster < ',num2str(th_numpix),' pix']})
            title('Isotropic Diffusion : Performance Dependence on Cluster Sizes')
            xlabel('AUCm - AUCi')
            ylabel('Counts')

            text(-0.4, 0.85*max([h1,h2]), ['\color{red}{n=',num2str(n),'}']);
            text(-0.4, 0.9*max([h1,h2]), ['\color{blue}{n=',num2str(numel(ind)-n),'}']);
        end


        % Modularity SKH Distribution split by large and small cluster pairs
        if(0)
            SK_AUCdiff = SK_101_5_lrg.kur_ROC_AUC_GS_accum - SK_101_5_lrg.im_ROC_AUC_1D_accum;

            th_numpix = 1000;
            ind = (SK_101_5_lrg.ClusterSizes(1,:) < th_numpix & SK_101_5_lrg.ClusterSizes(2,:) < th_numpix);
            n = numel(find(ind>0));
            h1 = hist(SK_AUCdiff(ind),bins);
            h2 = hist(SK_AUCdiff(~ind),bins);

            figure, hold on,
            bar(bins,h2,0.50,'b')
            bar(bins,h1,0.25,'r')
            xlim([-0.55,0.55])
            legend({['Both Clusters > ',num2str(th_numpix),' pix'],['Either Cluster < ',num2str(th_numpix),' pix']})
            title('Topographic Modularity : Performance Dependence on Cluster Sizes')
            xlabel('AUCm - AUCi')
            ylabel('Counts')

            text(-0.4, 0.85*max([h1,h2]), ['\color{red}{n=',num2str(n),'}']);
            text(-0.4, 0.9*max([h1,h2]), ['\color{blue}{n=',num2str(numel(ind)-n),'}']);
        end


        % Directly compare Modularity SKH to Isotropic Diffusion
        disp('Hmm, suspiciously, these (IsoDiff & ModSKH) distibutions look exactly the same.')
        disp('Compare SK & Iso more directly.')

        h = hist(SK_101_5_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_mid.kur_ROC_AUC_GS_accum, bins);

        [bins(1:10);bins(21:-1:12)]
        [h(1:10);h(21:-1:12)]


        disp('# Cluster Pairs where Mod SKH outperforms Isotropic Diffusion:')
        sum(h(12:21))

        disp('# Cluster Pairs where Isotropic Diffusion outperforms Mod SKH:')
        sum(h(1:10))

        disp('# Cluster Pairs where Isotropic Diffusion & Mod SKH perform same:')
        sum(h(11))




        % Directly compare Graph Laplacian to Isotropic Diffusion
        disp('Compare GL & Iso directly.')

        h = hist(GL_101_1_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_mid.kur_ROC_AUC_GS_accum, bins);

        [bins(1:10);bins(21:-1:12)]
        [h(1:10);h(21:-1:12)]


        disp('# Cluster Pairs where GL outperforms Isotropic Diffusion:')
        sum(h(12:21))

        disp('# Cluster Pairs where Isotropic Diffusion outperforms GL:')
        sum(h(1:10))

        disp('# Cluster Pairs where Isotropic Diffusion & GL perform same:')
        sum(h(11))




        % Directly compare Modularity NG to Isotropic Diffusion
        disp('Compare NG & Iso directly.')

        h = hist(NG_101_10_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_mid.kur_ROC_AUC_GS_accum, bins);

        [bins(1:10);bins(21:-1:12)]
        [h(1:10);h(21:-1:12)]


        disp('# Cluster Pairs where NG outperforms Isotropic Diffusion:')
        sum(h(12:21))

        disp('# Cluster Pairs where Isotropic Diffusion outperforms NG:')
        sum(h(1:10))

        disp('# Cluster Pairs where Isotropic Diffusion & NG perform same:')
        sum(h(11))



        % Directly compare Average Association to Isotropic Diffusion
        disp('Compare AA & Iso directly.')

        h = hist(AA_101_1_lrg.kur_ROC_AUC_GS_accum - Iso_101_3_mid.kur_ROC_AUC_GS_accum, bins);

        [bins(1:10);bins(21:-1:12)]
        [h(1:10);h(21:-1:12)]


        disp('# Cluster Pairs where AA outperforms Isotropic Diffusion:')
        sum(h(12:21))

        disp('# Cluster Pairs where Isotropic Diffusion outperforms AA:')
        sum(h(1:10))

        disp('# Cluster Pairs where Isotropic Diffusion & AA perform same:')
        sum(h(11))


























        %% Display Example Images that fall in certain bins in the Isotropic Diffusion AUCm-AUCi distribution

        th_numpix = 1000;
        ind1 = (Iso_101_3_mid.ClusterSizes(1,:) > th_numpix & Iso_101_3_mid.ClusterSizes(2,:) > th_numpix);
        ind2 = (SK_101_5_lrg.ClusterSizes(1,:) > th_numpix & SK_101_5_lrg.ClusterSizes(2,:) > th_numpix);
        ind3 = (NG_101_10_lrg.ClusterSizes(1,:) > th_numpix & NG_101_10_lrg.ClusterSizes(2,:) > th_numpix);
        ind4 = (GL_101_1_lrg.ClusterSizes(1,:) > th_numpix & GL_101_1_lrg.ClusterSizes(2,:) > th_numpix);
        ind5 = (AA_101_1_lrg.ClusterSizes(1,:) > th_numpix & AA_101_1_lrg.ClusterSizes(2,:) > th_numpix);
        %
        if ~isempty(find(ind1 ~= ind2))
            disp('Something doesnt match with number of and size of cluster pairs in Iso vs. SK.')
        end
        if ~isempty(find(ind1 ~= ind3))
            disp('Something doesnt match with number of and size of cluster pairs in Iso vs. NG.')
        end
        if ~isempty(find(ind1 ~= ind4))
            disp('Something doesnt match with number of and size of cluster pairs in Iso vs. GL.')
        end
        if ~isempty(find(ind1 ~= ind5))
            disp('Something doesnt match with number of and size of cluster pairs in Iso vs. AA.')
        end
        %
        Iso_best_filt = Iso_101_3_mid.kur_ROC_AUC_GS_accum(ind1);
        SK_best_filt = SK_101_5_lrg.kur_ROC_AUC_GS_accum(ind2);
        NG_best_filt = NG_101_10_lrg.kur_ROC_AUC_GS_accum(ind3);
        GL_best_filt = GL_101_1_lrg.kur_ROC_AUC_GS_accum(ind4);
        AA_best_filt = AA_101_1_lrg.kur_ROC_AUC_GS_accum(ind5);
        %
        Pix_Iso_filt = Iso_101_3_mid.im_ROC_AUC_1D_accum(ind1);
        Pix_SK_filt = SK_101_5_lrg.im_ROC_AUC_1D_accum(ind2);
        Pix_NG_filt = NG_101_10_lrg.im_ROC_AUC_1D_accum(ind3);
        Pix_GL_filt = GL_101_1_lrg.im_ROC_AUC_1D_accum(ind4);
        Pix_AA_filt = AA_101_1_lrg.im_ROC_AUC_1D_accum(ind5);
        %
        ClusterSize_Iso = Iso_101_3_mid.ClusterSizes(:,ind1);
        ClusterSize_SK =  SK_101_5_lrg.ClusterSizes(:,ind2);
        ClusterSize_NG =  NG_101_10_lrg.ClusterSizes(:,ind3);
        ClusterSize_GL =  GL_101_1_lrg.ClusterSizes(:,ind4);
        ClusterSize_AA =  AA_101_1_lrg.ClusterSizes(:,ind5);
        %
        ImgPatchID_Iso = Iso_101_3_mid.ImgPatchID(ind1);
        ImgPatchID_SK = SK_101_5_lrg.ImgPatchID(ind2);
        ImgPatchID_NG = NG_101_10_lrg.ImgPatchID(ind3);
        ImgPatchID_GL = GL_101_1_lrg.ImgPatchID(ind4);
        ImgPatchID_AA = AA_101_1_lrg.ImgPatchID(ind5);
        %
        GndTruthID_Iso = Iso_101_3_mid.GndTruthID(ind1);
        GndTruthID_SK = SK_101_5_lrg.GndTruthID(ind2);
        GndTruthID_NG = NG_101_10_lrg.GndTruthID(ind3);
        GndTruthID_GL = GL_101_1_lrg.GndTruthID(ind4);
        GndTruthID_AA = AA_101_1_lrg.GndTruthID(ind5);

        %
        ClusterPairID_Iso = Iso_101_3_mid.ClusterPairID(:,ind1);
        ClusterPairID_SK = SK_101_5_lrg.ClusterPairID(:,ind2);
        ClusterPairID_NG = NG_101_10_lrg.ClusterPairID(:,ind3);
        ClusterPairID_GL = GL_101_1_lrg.ClusterPairID(:,ind4);
        ClusterPairID_AA = AA_101_1_lrg.ClusterPairID(:,ind5);


    %     % find image patch to visualize where SK beats Iso.
    %     x2 = SK_best_filt - Iso_best_filt;
    %     [y2,i2] = sort(x2,'descend');


        % find image patch to visualize where Iso performs best.
        [y2,i2] = sort(Iso_best_filt,'descend');

        % (maybe) find image patch to visualize where Iso beats SK

        % (maybe) find image patch to visualize where SK & Iso perform about the same


        if(1)
            disp('These should all be the same')
            num = 10;
            Iso_best_filt(i2(1:num))
            SK_best_filt(i2(1:num))
            NG_best_filt(i2(1:num))
            GL_best_filt(i2(1:num))
            AA_best_filt(i2(1:num))
            %
            Pix_Iso_filt(i2(1:num))
            Pix_SK_filt(i2(1:num))
            Pix_NG_filt(i2(1:num))
            Pix_GL_filt(i2(1:num))
            Pix_AA_filt(i2(1:num))
            %
            ClusterSize_Iso(:,i2(1:num))
            ClusterSize_SK(:,i2(1:num))
            ClusterSize_NG(:,i2(1:num))
            ClusterSize_GL(:,i2(1:num))
            ClusterSize_AA(:,i2(1:num))
            %
            ImgPatchID_Iso(i2(1:num))
            ImgPatchID_SK(i2(1:num))
            ImgPatchID_NG(i2(1:num))
            ImgPatchID_GL(i2(1:num))
            ImgPatchID_AA(i2(1:num))
            %
            GndTruthID_Iso(i2(1:num))
            GndTruthID_SK(i2(1:num))
            GndTruthID_NG(i2(1:num))
            GndTruthID_GL(i2(1:num))
            GndTruthID_AA(i2(1:num))
            %
            ClusterPairID_Iso(:,i2(1:num))
            ClusterPairID_SK(:,i2(1:num))
            ClusterPairID_NG(:,i2(1:num))
            ClusterPairID_GL(:,i2(1:num))
            ClusterPairID_AA(:,i2(1:num))
        end

    end
    
    
    
    %% Visualize image patch where SK outperforms Iso by the most.
    
    
    if(1)
        
        
        dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/methodsComparePlots/'];
        if ~exist(dirImgSave,'dir')
            mkdir(dirImgSave);
        end
        
        
        % dirMatSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/']
        
        

        %ptch = '118015_ptch2'; % Mod_SKH performs much better than IsoDiff
        %ptch = '117054_ptch2'; % Iso Diff performs well here
        ptch = '100080_ptch2'; % Iso Diff performs well here too

        IsoK = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/IsoDiff/KurMC_',ptch,'_rM3_NF_60_0_ksmid.mat']);
        IsoE = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/IsoDiff/Evecs_',ptch,'_rM3.mat']);
        %
        SK_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/Mod_SKHAdj/KurMC_',ptch,'_rM5_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %SK_E =
        %
        NG_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/Mod_N&G/KurMC_',ptch,'_rM10_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %NG_E =
        %
        AA_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/AAnrm/KurMC_',ptch,'_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %AA_E =
        %
        GL_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/GLnrm/KurMC_',ptch,'_rM1_sDInf_sP0p2_NF_60_0_kslrg.mat']);
        %GL_E =


        GTFile = load([dirPre,'images/BSDS_patch/101x101_ds1/',ptch,'.mat']);
       

        if(1)
        
            % Instead of using each ground truth boundary individually, build
            % up one ground truth by the consensus of the 5 individuals
            bD_consensus = (GTFile.bD>1);
            bD2c = logical( blur( single(bD_consensus),1,'binom2' ) );
            bD3c = logical( blur( single(bD_consensus),1,'binom3' ) );
            bD4c = logical( blur( single(bD_consensus),1,'binom4' ) );
            bD5c = logical( blur( single(bD_consensus),1,'binom5' ) ); % note: 'power' increases with increased blurring


            bDc_blurs(:,:,1) = bD_consensus;
            bDc_blurs(:,:,2) = bD2c;
            bDc_blurs(:,:,3) = bD3c;
            bDc_blurs(:,:,4) = bD4c;
            bDc_blurs(:,:,5) = bD5c;

            
            % Number of cluster boundaries to cross image is likely
            % proportional to the amount of phase jump I can expect at the
            % boundaries since everyhting has to fit into 2pi.
            disp('As I cross the ground truth edges image horizontally or vertically,')
            disp(' the average number of boundaries I cross is ...')
            mean(sum(bD_consensus,2))
            mean(sum(bD_consensus,1))
            



            % As a sanity check, create a segmentation that consists of a
            % blurred version of the ground truth.  What do M & S look like?
            if(0)
                disp('Sanity Check: Blur 1')
                X = single(bD_consensus);
                [F,M.b1,S.b1] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                disp('Sanity Check: Blur 2')
                X = single(bD2c);
                [F,M.b2,S.b2] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                disp('Sanity Check: Blur 3')
                X = single(bD3c);
                [F,M.b3,S.b3] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                disp('Sanity Check: Blur 4')
                X = single(bD4c);
                [F,M.b4,S.b4] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                disp('Sanity Check: Blur 5clr')
                X = single(bD5c);
                [F,M.b5,S.b5] = compute_BoundaryGradientMetric(X,bDc_blurs);
                
                % Plot ratio of mean gradients for different methods (on boundary vs off)
                gg = [M.b1(:,1)./S.b1(:,1), ...
                      M.b2(:,1)./S.b2(:,1), ...
                      M.b3(:,1)./S.b3(:,1), ...
                      M.b4(:,1)./S.b4(:,1), ...
                      M.b5(:,1)./S.b5(:,1)];

                figure, 
                subplot(121), plot(gg,'LineWidth',2)
                title('Ratio of Mean Gradient Value On vs Off Boundaries')
                xlabel('blur')
                %legend({'pix','iso','sk','ng','gl','aa'}) 



                gg = [M.b1(:,1), ...
                      M.b2(:,1), ...
                      M.b3(:,1), ...
                      M.b4(:,1), ...
                      M.b5(:,1)];

                subplot(122), plot(gg,'LineWidth',2)
                title('Mean Gradient Value On Boundaries')
                xlabel('blur')
                legend({'b1','2','3','4','5'}) 
                
            end

            clear bD_consensus bD2c bD3c bD4c bD5c

            
            % PLOT #1: Visualize Gradients (phase,contrast or otherwise)
            % for different competitor methods and plot avg gradient/pixel
            % on boundaries vs. off for different boundary blurs.
            if(1)
                
                
                
                % Compute Mean & Std of Gradient Fields On & Off Ground Truth Boundaries
                X = GTFile.im;
                [F.im,M.im,S.im] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                X = visKurPhase_inHSV(IsoK.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end),IsoK.netParams.Ndims));
                [F.iso,M.iso,S.iso] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                X = visKurPhase_inHSV(NG_K.netParams.im, reshape(NG_K.metaCluster.phaseAtClk(:,end),NG_K.netParams.Ndims));
                [F.ng,M.ng,S.ng] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                X = visKurPhase_inHSV(SK_K.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end),SK_K.netParams.Ndims));
                [F.sk,M.sk,S.sk] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                X = visKurPhase_inHSV(GL_K.netParams.im, reshape(GL_K.metaCluster.phaseAtClk(:,end),GL_K.netParams.Ndims));
                [F.gl,M.gl,S.gl] = compute_BoundaryGradientMetric(X,bDc_blurs);
                %
                X = visKurPhase_inHSV(AA_K.netParams.im, reshape(AA_K.metaCluster.phaseAtClk(:,end),AA_K.netParams.Ndims));
                [F.aa,M.aa,S.aa] = compute_BoundaryGradientMetric(X,bDc_blurs);
                

                % Find max value for y-axis of error bar plots (M_mn+M_std or S_mn+S_std)
                del_mn_std = [sum(M.im,2); sum(M.iso,2)./pi; sum(M.ng,2)./pi; sum(M.sk,2)./pi; sum(M.gl,2)./pi; sum(M.aa,2)./pi;...
                    sum(S.im,2); sum(S.iso,2)./pi; sum(S.ng,2)./pi; sum(S.sk,2)./pi; sum(S.gl,2)./pi; sum(S.aa,2)./pi];


                % Figure with Subplots !
                Hc=figure;
                %
                % Ground Truth Boundaries Drawn by Different Human Subjects
                figure(Hc), subplot(441), imagesc(GTFile.bD), colormap('jet'), colorbar, axis square, title('Indiv. Gnd Truth','FontSize',18,'FontWeight','Bold')
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze
                %
                % Consensus Ground Truth Boundaries (where 2 or more human subjects agree)
                figure(Hc), subplot(442), imagesc(GTFile.bD>1), colormap('bone'), colorbar, axis square, title('Consensus gT','FontSize',18,'FontWeight','Bold')
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze
                %
                % Model 1: Just Image Pixels
                figure(Hc), subplot(443), imagesc(F.im), colormap('jet'), cb=colorbar; set(get(cb,'xlabel'),'string','\Delta\Theta')
                axis square, title('Image Pixels','FontSize',18,'FontWeight','Bold'),
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze
                subplot(447), hold on, errorbar(M.im(:,1),M.im(:,2),'b.-'), errorbar(S.im(:,1),S.im(:,2),'r.-')
                set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 size(bDc_blurs,3)+0.5 0 1.05*max(del_mn_std)])
                %
                % Model 2: Isotropic Diffusion (Oscillator Relaxation)
                figure(Hc), subplot(444), imagesc(F.iso), colormap('jet'), colorbar, axis square, title('Isotropic Diffusion','FontSize',18,'FontWeight','Bold') 
                set(gca,'XTick',[],'YTick',[]), 
                subplot(448), hold on, errorbar(M.iso(:,1)./pi,M.iso(:,2)./pi,'b.-'), errorbar(S.iso(:,1)./pi,S.iso(:,2)./pi,'r.-')
                set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 size(bDc_blurs,3)+0.5 0 1.05*max(del_mn_std)])
                legend({'on boundary','off boundary'})
                %
                % ASIDE: Plot Image Patch just for Observer Orientation.
                subplot(4,4,[5,6]), imagesc(GTFile.im), colormap('bone'), axis square,
                set(gca,'XTick',[],'YTick',[]), freezeColors,
                ptch_str = ptch;
                ptch_str(ptch_str=='_')=' ';
                title(['Im Patch : ',ptch_str],'FontSize',18,'FontWeight','Bold') 
                %
                % Model 3: Newman & Girvan Modularity (Oscillator Relaxation)
                figure(Hc), subplot(449), imagesc(F.ng), colormap('jet'), colorbar, axis square, title('Modularity NG','FontSize',18,'FontWeight','Bold')
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze 
                subplot(4,4,13), hold on, errorbar(M.ng(:,1)./pi,M.ng(:,2)./pi,'b.-'), errorbar(S.ng(:,1)./pi,S.ng(:,2)./pi,'r.-')
                set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 size(bDc_blurs,3)+0.5 0 1.05*max(del_mn_std)])
                xlabel('gT boundary blur')
                ylabel({'mean \Delta','(% of D.R.)'})
                %
                % Model 4: SKH Topographic Modularity (Oscillator Relaxation)
                figure(Hc), subplot(4,4,10), imagesc(F.sk), colormap('jet'), colorbar, axis square, title('Modularity SK','FontSize',18,'FontWeight','Bold') 
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze
                subplot(4,4,14), hold on, errorbar(M.sk(:,1)./pi,M.sk(:,2)./pi,'b.-'), errorbar(S.sk(:,1)./pi,S.sk(:,2)./pi,'r.-')
                set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 size(bDc_blurs,3)+0.5 0 1.05*max(del_mn_std)])
                %
                % Model 5: Graph Laplacian (Oscillator Relaxation)
                figure(Hc), subplot(4,4,11), imagesc(F.gl), colormap('jet'), colorbar, axis square, title('Graph Laplacian','FontSize',18,'FontWeight','Bold') 
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze
                subplot(4,4,15), hold on, errorbar(M.gl(:,1)./pi,M.gl(:,2)./pi,'b.-'), errorbar(S.gl(:,1)./pi,S.gl(:,2)./pi,'r.-')
                set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 size(bDc_blurs,3)+0.5 0 1.05*max(del_mn_std)])
                %
                % Model 6: Average Association (Oscillator Relaxation)
                figure(Hc), subplot(4,4,12), imagesc(F.aa), colormap('jet'), colorbar, axis square, title('Average Association','FontSize',18,'FontWeight','Bold') 
                set(gca,'XTick',[],'YTick',[]), freezeColors, cbfreeze
                subplot(4,4,16), hold on, errorbar(M.aa(:,1)./pi,M.aa(:,2)./pi,'b.-'), errorbar(S.aa(:,1)./pi,S.aa(:,2)./pi,'r.-')
                set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 size(bDc_blurs,3)+0.5 0 1.05*max(del_mn_std)])
                
                
                saveGoodImg(Hc,[dirImgSave,ptch,'_visGradients'],sizeGoodIm)
                close(Hc)

            end
            
            
            
            
            % PLOT #2: Plot 
            if(1)
            
                % Plot ratio of mean gradients for different methods (on boundary vs off)
                gg1 = [M.im(:,1)./S.im(:,1), ...
                      M.iso(:,1)./S.iso(:,1), ...
                      M.sk(:,1)./S.sk(:,1), ...
                      M.ng(:,1)./S.ng(:,1), ...
                      M.gl(:,1)./S.gl(:,1), ...
                      M.aa(:,1)./S.aa(:,1)];

                Hl = figure; 
                subplot(121), plot(gg1,'LineWidth',2)
                title('Ratio of Mean Gradient Value On vs Off Boundaries','Fontweight','Bold','FontSize',20)
                xlabel('blur','Fontweight','Bold','FontSize',18)
                ylabel('M_{B}:/:M_{~B}','Fontweight','Bold','FontSize',18)
                set(gca,'FontSize',16,'FontWeight','Bold')
                %legend({'pix','iso','sk','ng','gl','aa'}) 
                text(size(gg1,1)-1, max(gg1(:)), [ptch_str],'Fontweight','Bold','FontSize',16,'HorizontalAlignment','Center')



                gg2 = [M.im(:,1), ...
                      M.iso(:,1)./pi, ...
                      M.sk(:,1)./pi, ...
                      M.ng(:,1)./pi, ...
                      M.gl(:,1)./pi, ...
                      M.aa(:,1)./pi];

                subplot(122), plot(gg2,'LineWidth',2)
                title('Mean Gradient Value On Boundaries','Fontweight','Bold','FontSize',20)
                xlabel('blur','Fontweight','Bold','FontSize',18)
                ylabel('M_{B} (% of Dynamic Range)','Fontweight','Bold','FontSize',18)
                set(gca,'FontSize',16,'FontWeight','Bold')
                legend({'pix','iso','sk','ng','gl','aa'}) 
                
                saveGoodImg(Hl,[dirImgSave,ptch,'_grads_OnNoffBoundary'],sizeGoodIm)
                close(Hl)

            end
            
            
            
        end
        
        
        
        
        
        
        
        
        
        
        
        
        % find the entries indexed by i2 that come from this image patch.0        pp = strfind(ImgPatchID_Iso(i2),ptch);
        ptch_ind = [];
        for r = 1:numel(pp)
            if ~isempty(pp{r})
                ptch_ind = [ptch_ind,r];
            end
        end

        ImgPatchID_Iso(i2(ptch_ind))



        % further sort data from this image patch by the ground truth it is associated with.
        [l,gt_ind] = sort(GndTruthID_Iso(i2(ptch_ind)))








        gg = [Pix_Iso_filt(i2(ptch_ind(gt_ind))); ...
              Iso_best_filt(i2(ptch_ind(gt_ind))); ...
               SK_best_filt(i2(ptch_ind(gt_ind))); ...
               NG_best_filt(i2(ptch_ind(gt_ind))); ...
               GL_best_filt(i2(ptch_ind(gt_ind))); ...
               AA_best_filt(i2(ptch_ind(gt_ind)))]


        if(0)

            % Plot AUC distributions for different methods as well as final phase.
            figure,
            subplot(4,4,[1,2,5,6]), hold on, 
            % replotting distributions of (AUCm-AUCi) distributions for isotropic diffusion and 4 anisotropic diffusion models.
        %     plot(bins, normlze(hIso_101_3_mid),'g','LineWidth',2)       % IsoDiff Kuramoto
        %     plot(bins, normlze(hSK_101_5_lrg),'r','LineWidth',2)         % SK Kuramoto
        %     plot(bins, normlze(hNG_101_10_lrg),'c','LineWidth',2)       % NG Kuramoto
        %     plot(bins, normlze(hGL_101_1_lrg),'m','LineWidth',2)         % GL Kuramoto
        %     %plot(bins, normlze(hGLev_101_3_lrg),'m--','LineWidth',2)  % GL 2nd eigenvector
        %     plot(bins, normlze(hAA_101_1_lrg),'y','LineWidth',2)         % AA Kuramoto

            plot(repmat(bins',1,6),...
                       [zeros(numel(bins),1),     ...                              % Just Image Pixels
                        normlze(hIso_101_3_mid)', ...                              % IsoDiff Kuramoto
                        normlze(hSK_101_5_lrg)',  ...                              % SK Kuramoto
                        normlze(hNG_101_10_lrg)', ...                              % NG Kuramoto
                        normlze(hGL_101_1_lrg)', ...                               % GL Kuramoto
                        normlze(hAA_101_1_lrg)'],'LineWidth',2)                    % AA Kuramoto
            %plot(bins, normlze(hGLev_101_3_lrg),'m--','LineWidth',2)              % GL 2nd eigenvector
            %
            plot([0 0], [0 0.4], 'k--','LineWidth',1.1)
            xlabel('AUC_{method} - AUC_{img}','FontSize',18,'FontWeight','Bold')
            ylabel('counts (segment pairs)','FontSize',18,'FontWeight','Bold')
            title('Optimized Performance Distributions','FontSize',16,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            legend({'pix','iso','sk','ng','gl','aa'})  % 'GL Ev2'
            xlim([-0.6 0.6])
            ylim([0 0.5])
            %
            herrorbar(mIso_101_3_mid, 0.47, sIso_101_3_mid, 'g.')
            herrorbar(mSK_101_5_lrg, 0.46, sSK_101_5_lrg, 'r.')
            herrorbar(mNG_101_10_lrg, 0.45, sNG_101_10_lrg, 'c.')
            herrorbar(mGL_101_3_lrg, 0.44, sGL_101_3_lrg, 'm.')
            herrorbar(mAA_101_1_lrg, 0.43, sAA_101_1_lrg, 'y.')
            text(0.1,0.5,'Mean & STD','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')


            subplot(4,4,[9,10]), hold on, 
            plot(gg' - repmat(gg(:,1)',1,6),'.','LineWidth',2)
            legend({'pix','iso','sk','ng','gl','aa'})   
            ylim([-0.6 0.6])
            ylabel('AUC')
            xlabel('Cluster Pair Index')
            title(ptch)
            %
            errorbar( size(gg,2)+2, mean(Pix_Iso_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(Pix_Iso_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'bd', 'LineWidth',2 )
            errorbar( size(gg,2)+4, mean(Iso_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(Iso_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'gd', 'LineWidth',2 )
            errorbar( size(gg,2)+6, mean(SK_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(SK_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'rd', 'LineWidth',2 )
            errorbar( size(gg,2)+8, mean(NG_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(NG_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'cd', 'LineWidth',2 )
            errorbar( size(gg,2)+10, mean(GL_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(GL_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'md', 'LineWidth',2 )
            errorbar( size(gg,2)+12, mean(AA_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(AA_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'yd', 'LineWidth',2 )

            view([90 90]);

            % (1). Image Pixels Strawman
            subplot(4,4,3), imagesc(IsoE.netParams.im), colormap('bone'), axis square, title('Im Pix')
            set(gca,'XTick',[],'YTick',[])

            % (2). Isotropic Diffusion
            subplot(4,4,4), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end), IsoK.netParams.Ndims) ) ), colormap('hsv'), axis square, title('Iso Diff.')
            set(gca,'XTick',[],'YTick',[])

            % (3). Modularity SKH
            subplot(4,4,7), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end), SK_K.netParams.Ndims) ) ), colormap('hsv'), axis square, title('Mod SK.')
            set(gca,'XTick',[],'YTick',[])

            % (4). Graph Laplacian
            subplot(4,4,8), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(GL_K.metaCluster.phaseAtClk(:,end), GL_K.netParams.Ndims) ) ), colormap('hsv'), axis square, title('GL')
            set(gca,'XTick',[],'YTick',[])

            % (5). Modularity N&G
            subplot(4,4,11), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(NG_K.metaCluster.phaseAtClk(:,end), NG_K.netParams.Ndims) ) ), colormap('hsv'), axis square, title('Mod N&G')
            set(gca,'XTick',[],'YTick',[])

            % (6). Average Association
            subplot(4,4,12), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(AA_K.metaCluster.phaseAtClk(:,end), AA_K.netParams.Ndims) ) ), colormap('hsv'), axis square, title('AA')
            set(gca,'XTick',[],'YTick',[])

            % Plot ground truths
            for kk = 1:numel(IsoK.netParams.gT)
                subplot(4, numel(IsoK.netParams.gT), 3*numel(IsoK.netParams.gT)+kk), imagesc(IsoK.netParams.gT{kk}), axis square
                set(gca,'XTick',[],'YTick',[])
                xlabel(['gT#',num2str(kk)])
            end

        end

        % % % % % %



        % Can I visualize which 2 segments in a ground truth a set of 6 points is looking at?
        cmap_rwb = rd_plotColorbar('redwhiteblue',128);

        figure
        subplot(4,1,1), hold on, grid on
        plot(gg' - repmat(gg(:,1)',1,6) ,'.','MarkerSize',30) % 
        legend({'pix','iso','sk','ng','gl','aa'})   
        ylim([0 0.5])
        ylabel('AUCm-AUCi','FontSize',16,'FontWeight','Bold')
        xlabel('Cluster Pair Index','FontSize',16,'FontWeight','Bold')
        set(gca,'XTick',1:(size(gg,2)+6),'XTickLabel',[1:size(gg,2),zeros(1,6)],'FontSize',14,'FontWeight','Bold')
        title(ptch)

        errorbar( size(gg,2)+1, mean(Pix_Iso_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(Pix_Iso_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'b.', 'LineWidth',3,'MarkerSize',30 )
        errorbar( size(gg,2)+2, mean(Iso_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(Iso_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'g.', 'LineWidth',3,'MarkerSize',30 )
        errorbar( size(gg,2)+3, mean(SK_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(SK_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'r.', 'LineWidth',3,'MarkerSize',30 )
        errorbar( size(gg,2)+4, mean(NG_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(NG_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'c.', 'LineWidth',3,'MarkerSize',30 )
        errorbar( size(gg,2)+5, mean(GL_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(GL_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'm.', 'LineWidth',3,'MarkerSize',30 )
        errorbar( size(gg,2)+6, mean(AA_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), std(AA_best_filt(i2(ptch_ind))-Pix_Iso_filt(i2(ptch_ind))), 'y.', 'LineWidth',3,'MarkerSize',30 )



        for kk = 1:numel(gt_ind)

            GT = int16(IsoK.netParams.gT{GndTruthID_Iso(i2(ptch_ind(gt_ind(kk))))});
            GT(GT==ClusterPairID_Iso(1,i2(ptch_ind(gt_ind(kk))))) = 1000;
            GT(GT==ClusterPairID_Iso(2,i2(ptch_ind(gt_ind(kk))))) = -1000;
            GT(abs(GT)<1000) = 0;

            subplot(4,numel(gt_ind),numel(gt_ind)+kk)
            imagesc(GT), colormap(cmap_rwb), caxis([-1000 1000]), freezeColors
            axis square
            title(num2str(kk),'FontSize',16,'FontWeight','Bold')
            set(gca,'XTick',[],'YTick',[],'FontSize',14,'FontWeight','Bold')
            if(kk==1)
                ylabel('gT#','FontSize',16,'FontWeight','Bold')
            end
            xlabel(['\color{blue}{',num2str(Pix_Iso_filt(i2(ptch_ind(kk))),2),'}'])

        end


        % (1). Image Pixels Strawman
        subplot(4,6,13), imagesc(IsoE.netParams.im), colormap('bone'), caxis([0 1]), axis square, 
        title('\color{blue}{Im Pix}','FontSize',16,'FontWeight','Bold'), freezeColors
        set(gca,'XTick',[],'YTick',[])

        % (2). Isotropic Diffusion
        subplot(4,6,14), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end), IsoK.netParams.Ndims) ) ),
        colormap('hsv'), caxis([0 2*pi]), axis square, 
        title('\color{green}{Iso Diff}','FontSize',16,'FontWeight','Bold'), freezeColors
        set(gca,'XTick',[],'YTick',[])

        % (3). Modularity SKH
        subplot(4,6,15), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end), SK_K.netParams.Ndims) ) ), 
        colormap('hsv'), caxis([0 2*pi]), axis square, 
        title('\color{red}{Mod SK}','FontSize',16,'FontWeight','Bold'), freezeColors
        set(gca,'XTick',[],'YTick',[])

        % (4). Modularity N&G
        subplot(4,6,16), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(NG_K.metaCluster.phaseAtClk(:,end), NG_K.netParams.Ndims) ) ), 
        colormap('hsv'), caxis([0 2*pi]), axis square, 
        title('\color{cyan}{Mod N&G}','FontSize',16,'FontWeight','Bold'), freezeColors
        set(gca,'XTick',[],'YTick',[])

        % (5). Graph Laplacian
        subplot(4,6,17), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(GL_K.metaCluster.phaseAtClk(:,end), GL_K.netParams.Ndims) ) ), 
        colormap('hsv'), caxis([0 2*pi]), axis square, 
        title('\color{magenta}{GL}','FontSize',16,'FontWeight','Bold'), freezeColors
        set(gca,'XTick',[],'YTick',[])

        % (6). Average Association
        subplot(4,6,18), imagesc( visKurPhase_inHSV( IsoE.netParams.im, reshape(AA_K.metaCluster.phaseAtClk(:,end), AA_K.netParams.Ndims) ) ), 
        colormap('hsv'), caxis([0 2*pi]), axis square, 
        title('\color{yellow}{AA}','FontSize',16,'FontWeight','Bold'), freezeColors
        set(gca,'XTick',[],'YTick',[])


        % Plot Histograms of the Pixel and Phase Distributions.
        colorsHSV = colormap('hsv');
        colorsHSV = colorsHSV(1:2:end,:);
        colorsBone = colormap('bone');
        colorsBone = colorsBone(1:2:end,:);
        %
        histBinsCirc = linspace(0,2*pi,64);
        histBinsLinr = linspace(0,1,64);
        %
        % for image pixels
        xy = IsoE.netParams.im ;
        [n,x] = hist(xy(:),histBinsLinr);
        subplot(4,6,19)
        hBar = bar(x,n,'hist');
        set(hBar,'FaceVertexCData',colorsBone);
        xlim([0 1])
        %
        % for isotropic diffusion
        xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end), SK_K.netParams.Ndims) ) ;
        [n,x] = hist(xy(:),histBinsCirc);
        subplot(4,6,20)
        hBar = bar(x,n,'hist');
        set(hBar,'FaceVertexCData',colorsHSV);
        xlim([0 2*pi])
        %
        % for modularity skh
        xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end), SK_K.netParams.Ndims) ) ;
        [n,x] = hist(xy(:),histBinsCirc);
        subplot(4,6,21)
        hBar = bar(x,n,'hist');
        set(hBar,'FaceVertexCData',colorsHSV);
        xlim([0 2*pi])
        %
        % for modularity n&g
        xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(NG_K.metaCluster.phaseAtClk(:,end), NG_K.netParams.Ndims) ) ;
        [n,x] = hist(xy(:),histBinsCirc);
        subplot(4,6,22)
        hBar = bar(x,n,'hist');
        set(hBar,'FaceVertexCData',colorsHSV);
        xlim([0 2*pi])
        %
        % for graph laplacian
        xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(GL_K.metaCluster.phaseAtClk(:,end), GL_K.netParams.Ndims) ) ;
        [n,x] = hist(xy(:),histBinsCirc);
        subplot(4,6,23)
        hBar = bar(x,n,'hist');
        set(hBar,'FaceVertexCData',colorsHSV);
        xlim([0 2*pi])
        %
        % for average association
        xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(AA_K.metaCluster.phaseAtClk(:,end), AA_K.netParams.Ndims) ) ;
        [n,x] = hist(xy(:),histBinsCirc);
        subplot(4,6,24)
        hBar = bar(x,n,'hist');
        set(hBar,'FaceVertexCData',colorsHSV);
        xlim([0 2*pi])
        
        
        
        
        
        
        % Plot a 3D scatter plot with x & y position of oscillator and phase as z location
        if(0)
        X = repmat([1:SK_K.netParams.Ndims(1)],SK_K.netParams.Ndims(2),1);
            X1 = reshape(X,1,SK_K.netParams.N);
            X2 = reshape(X',1,SK_K.netParams.N);

            xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end), SK_K.netParams.Ndims) ) ;
            yz = reshape(xy,1,SK_K.netParams.N);
            ind = (yz < 1.0 | yz > 5.5);
            figure, colormap('hsv'), caxis([0 2*pi]), zlim([0 2*pi])
            scatter3([X1(ind),0,0], [X2(ind),101,101], [yz(ind),0,2*pi], 10, [yz(ind),0,2*pi]), title('SK'), axis ij

            xy = visKurPhase_inHSV( IsoE.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end), IsoK.netParams.Ndims) ) ;
            figure, colormap('hsv'), scatter3(X1, X2, reshape(xy,1,IsoK.netParams.N),10,reshape(xy,1,IsoK.netParams.N)), title('Iso'), axis ij

            xy = IsoE.netParams.im;
            figure, colormap('jet'), scatter3(X1, X2, reshape(xy,1,SK_K.netParams.N),10,reshape(xy,1,SK_K.netParams.N)), title('SK'), axis ij
        end
        

        ClusterSize_Iso(:,i2(ptch_ind(gt_ind)))
        ClusterSize_SK(:,i2(ptch_ind(gt_ind)))
        ClusterSize_NG(:,i2(ptch_ind(gt_ind)))
        ClusterSize_GL(:,i2(ptch_ind(gt_ind)))
        ClusterSize_AA(:,i2(ptch_ind(gt_ind)))
        %
        %
        GndTruthID_Iso(i2(ptch_ind(gt_ind)))
        GndTruthID_SK(i2(ptch_ind(gt_ind)))
        GndTruthID_NG(i2(ptch_ind(gt_ind)))
        GndTruthID_GL(i2(ptch_ind(gt_ind)))
        GndTruthID_AA(i2(ptch_ind(gt_ind)))
        %
        ClusterPairID_Iso(:,i2(ptch_ind(gt_ind)))
        ClusterPairID_SK(:,i2(ptch_ind(gt_ind)))
        ClusterPairID_NG(:,i2(ptch_ind(gt_ind)))
        ClusterPairID_GL(:,i2(ptch_ind(gt_ind)))
        ClusterPairID_AA(:,i2(ptch_ind(gt_ind)))



        ClusterSize_Iso(:,i2(ptch_ind(gt_ind)))
        GndTruthID_Iso(i2(ptch_ind(gt_ind)))
        ClusterPairID_Iso(:,i2(ptch_ind(gt_ind)))
    end
    
    
    
    
    
    
    
    
    
    keyboard
    
    

end % end main function

function xn = normlze(x)

    xn = x ./ sum(x);

end
