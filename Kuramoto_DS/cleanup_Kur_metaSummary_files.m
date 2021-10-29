function cleanup_Kur_metaSummary_files(start, finish)


fileType = 'BSDS_patch';
fileSize = '51x51_ds1';

[dirPre,sizeGoodIm] = onCluster;
dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];


Tscale = {'1','0p5','0p1'};


disp('Starting out with this many files')
filesKur = dir([dirKur,'Kur_metaSummary*']) % can loop through these later


% Loop through filesKur to get rid of repeats (i.e. only keep Kur_metaSummary file containing larges number of files)


% iids1 = imgList('all');

% Note there are only 490 of these.  And in the sort_BSDS_output_files function
iids = [100007 100039 100075 100080 100098 100099 10081 101027 101084 101085 101087 102061 102062 103006 103029 103041 103078 ...
    104010 104022 104055 105019 105025 105027 105053 106005 106020 106024 106025 106047 107014 107045 107072 108004 108005 108036 ...
    108041 108069 108070 108073 109034 117025 117054 118015 118020 118031 118035 118072 119082 120003 12003 120093 12074 12084 ...
    122048 123057 123074 124084 126007 126039 128035 130014 130026 130034 130066 134008 134035 134049 134052 134067 135037 135069 ...
    138032 138078 140006 140055 140075 140088 14037 14085 14092 141012 141048 143090 144067 145014 145053 145059 145079 145086 146074 ...
    147021 147062 147077 147080 147091 148026 148089 15004 15011 15062 15088 151087 153077 153093 155060 156054 156065 156079 157032 ...
    157036 157055 157087 159002 159008 159022 159029 159045 159091 160006 16004 160067 160068 16052 16068 16077 161045 161062 163004 ...
    163014 163062 163085 163096 164046 164074 166081 167062 167083 168084 169012 170054 170057 17067 172032 173036 175032 175043 175083 ...
    176019 176035 176039 176051 178054 179084 181018 181021 181079 181091 182053 183055 183066 183087 185092 187003 187029 187039 187058 ...
    187071 187083 187099 188005 188025 188063 188091 189003 189006 189011 189013 189029 189080 189096 19021 196015 196027 196040 196062 ...
    196073 196088 197017 198004 198023 198054 198087 20008 20069 201080 2018 202000 202012 206062 206097 207038 207049 207056 208001 ...
    208078 209021 209070 2092 210088 21077 216041 216053 216066 216081 217013 217090 219090 220003 220075 22013 22090 22093 223004 ...
    223060 223061 225017 225022 226022 226033 226043 226060 227040 227046 227092 228076 229036 230063 230098 23025 23050 23080 23084 ...
    231015 232038 232076 235098 236017 236037 238011 238025 239007 239096 24004 24063 24077 241004 241048 242078 243095 245051 246009 ...
    246016 246053 247003 247012 247085 249021 249061 249087 250047 250087 25098 253016 253027 253036 253055 253092 254033 254054 257098 ...
    258089 259060 260058 260081 26031 267036 268002 268048 268074 27059 271008 271031 271035 274007 277053 277095 279005 28075 28083 28096 ...
    281017 285022 285036 285079 286092 288024 289011 290035 29030 291000 292066 293029 295087 296007 296028 296058 296059 299086 299091 ... 
    300091 301007 302003 302008 302022 304034 304074 306005 306051 306052 3063 309004 309040 3096 310007 311068 311081 314016 317043 317080 ...
    323016 326025 326038 326085 33039 33044 33066 334025 335088 335094 344010 346016 347031 35008 35010 35028 35049 35058 35070 35091 351093 ...
    353013 36046 361010 361084 365025 365072 365073 368016 368037 368078 370036 37073 372019 372047 374020 374067 376001 376020 376043 376086 ...
    38082 38092 384022 384089 385022 385028 385039 388006 388016 388018 388067 393035 41004 41006 41025 41029 41033 41069 41085 41096 42012 ...
    42044 42049 42078 43033 43051 43070 43074 43083 45000 45077 45096 46076 48017 48025 48055 49024 5096 51084 54005 54082 55067 55073 55075 ...
    56028 58060 59078 60079 6046 61034 61060 61086 62096 64061 65010 65019 65033 65074 65084 65132 66039 66053 66075 67079 68077 69000 69007 ...
    69015 69020 69022 69040 70011 70090 71046 71076 71099 76002 76053 77062 78004 78019 78098 79073 80085 80090 80099 8023 8049 8068 81066 ...
    81090 81095 8143 85048 86000 86016 86068 87015 87046 87065 89072 90076 92014 92059 94079 94095 95006 97010 97017 97033];


for i = start:finish % loop through all images in BSDS   
    
    disp(['File ',num2str(i),' / ',num2str(numel(iids)),'. Doing Between ',num2str(start),' & ',num2str(finish)])
    
    for j = 1:10 % number of patches extracted from each image
        
        j
        
        for k = 1:numel(Tscale)
    
            files = dir([dirKur,'Kur_metaSummary_BSDS_patch_',num2str(iids(i)),'_ptch',num2str(j),'__tscale',Tscale{k},'*'])
            
            % If there exist more than 1 Kur_metaSummary file, keep only the one containing the most files.
            for bj = 1:(numel(files)-1)
                delete([dirKur,files(bj).name])
            end
        
        end
    
    
    end
    
end


disp('Finishing with this many files')
filesKur = dir([dirKur,'Kur_metaSummary*']) % can loop through these later