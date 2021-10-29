function reduce_footprint_of_ModNG_mat_files(rM,ks,blur_flg)



[dirPre,sizeGoodIm] = onCluster;


if(blur_flg)
    blur_tag_M = '_blur_sig2';
else
    blur_tag_M = ''; % if we are not blurring.
end

specOrKur = 'spectral'; % 'Kur_PIF_Fourier1'

dirMatIn = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',specOrKur,'/Mod_N&G/'];

if strmatch('Kur_PIF_Fourier1',specOrKur)
    files = dir([dirMatIn,'KurMC_*_rM',rM,'_*_ks',ks,'.mat']);
elseif strmatch('spectral',specOrKur)
    files = dir([dirMatIn,'Evecs_*_rM',rM,'_*.mat']);
else
    disp(['Neither Spectral Nor Kur_PIF_Fourier1.  Check the specOrKur variable.'])
end

for i = 1:numel(files)
    
    disp([num2str(i),' / ',num2str(numel(files))])
    
    if files(i).bytes < 100e6
        
        disp('File less than 100MB. Ignoring it because its not taking up too much space.')
        [dirMatIn,files(i).name]
        
    else

        try

            load([dirMatIn,files(i).name])

            % remove two matrices which are taking up ALOT of harddrive space.
            netParams = rmfield(netParams,'Q');
            netParams = rmfield(netParams,'Wdist');
            
            
            % Save the mat file without those two large matrices.
            if strmatch('Kur_PIF_Fourier1',specOrKur)
                save([dirMatIn,files(i).name],'MC','metaCluster','kurParams','kurflags','netParams','netflags')
            elseif strmatch('spectral',specOrKur)
                save([dirMatIn,files(i).name],'MC','EVecsML','EValsML','EVecPM','netParams','netflags')
            else
                disp(['Neither Spectral Nor Kur_PIF_Fourier1.  Check the specOrKur variable.'])
            end
            

        catch

            disp('Must have already done this. Or File is Corrupt')
            [dirMatIn,files(i).name]

        end
    
    
    end

end