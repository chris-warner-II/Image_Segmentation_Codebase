
%
%
% 


methods = {'Q0'};
Ks = [100];
Ws = [0];

for m = 1:numel(methods)

    % Directory where mat files containing pairwise phase distance measures are located.
    dirDatIn = ['./output/explore_modularity_via_handcrafted_networks/pairwise_Kur_phase_diff/',methods{m},'/'];

    rpt_files = dir([dirDatIn,methods{m},'pairwise_Kur_phase_diff_Ks',num2str(Ks),'_sigW',num2str(Ws),'_rpt*.mat']);
    R = numel(rpt_files);
    
    
    % Initialize d_sum structure to be zero matrices of the right size.
    load([dirDatIn,rpt_files(1).name]);
    sz = size(d.mean_clust_pair.P1);
    num_tg = size(d.stats.traj_times1,1);   % number of time points or tp chunks grabbed from sim.
    num_clust = sz(1);                      % number of clusters embedded in network.
    num_tchunks = numel(d.stats.traj_tpts); % number of different chunkings of time
    %
    for p = 1:num_tchunks
        eval([ 'd_sum.mean_clust_pair.P',num2str(p),' = zeros([sz,R]);' ])
        eval([ 'd_sum.std_clust_pair.P',num2str(p),' = zeros([sz,R]);' ])
        %
        eval([ 'd_sum.mean_clust_pair.V',num2str(p),' = zeros([sz,R]);' ])
        eval([ 'd_sum.std_clust_pair.V',num2str(p),' = zeros([sz,R]);' ])
    end

    
    
    % Sum up various pairwise distance measures over a number of repeat_files made from multiple simulations.
    for r = 1:numel(rpt_files)
        load([dirDatIn,rpt_files(r).name])
        for p = 1:num_tchunks
            eval([ 'd_sum.mean_clust_pair.P',num2str(p),'(:,:,:,',num2str(r),') = d.mean_clust_pair.P',num2str(p),';' ])
            eval([ 'd_sum.std_clust_pair.P',num2str(p),'(:,:,:,',num2str(r),') = d.std_clust_pair.P',num2str(p),';' ])
            eval([ 'd_sum.mean_clust_pair.V',num2str(p),'(:,:,:,',num2str(r),') = d.mean_clust_pair.V',num2str(p),';' ])
            eval([ 'd_sum.std_clust_pair.V',num2str(p),'(:,:,:,',num2str(r),') = d.std_clust_pair.V',num2str(p),';' ])
        end
    end
    %
    d_sum.stats = d.stats;
    
    clear d
    
    
    squeeze(d_sum.mean_clust_pair.P1(num_clust,1,:,:)) % matrix of means one cluster pair across time points and repeats
    squeeze(d_sum.std_clust_pair.P1(num_clust,1,:,:))  % matrix of stds one cluster pair across time points and repeats
    
    
    
    %
    
    
    
    
    
    % HAVE TO TAKE MEANS AND STD'S MAX'S AND MIN'S OF PAIRWISE PHASE DIFFERENCES IN CLUSTERS (IN),
    % ACROSS CLUSTER BOUNDARIES (OUT), AND FOR ALL PAIRS (TOT) IN NETWORK...
    %
    
    
        
    for p = 1:num_tchunks % Loop through different time chunking methods.

        eval([' dmm = mean(d_sum.mean_clust_pair.P',num2str(p),',4); '])   % dmm is TEMP VAR ONLY.
        % *** dmm is MEAN across simulation repeats of MEAN pairwise distance between oscillators
        %     in cluster i and cluster j (size is num_clust x num_clust x num_time_pts)
        %
        eval([' dsm = std(d_sum.mean_clust_pair.P',num2str(p),',[],4); ']) % dsm is TEMP VAR ONLY.
        % *** dsm is STD across simulation repeats of MEAN pairwise distance between oscillators 
        %     in cluster i and cluster j (size is num_clust x num_clust x num_time_pts)
        %
        eval([' dms = mean(d_sum.std_clust_pair.P',num2str(p),',4); '])    % dms is TEMP VAR ONLY.
        % *** dms is MEAN across simulation repeats of STD of pairwise distance between oscillators
        %     in cluster i and cluster j (size is num_clust x num_clust x num_time_pts)
        
        
        for t = 1:num_tg % Loop through time points where time chunking happened and compute statistics AT EACH TIME POINT.
            %
            % % (1). MEAN & STD across clusters (in, out or tot network) of MEAN across simulation repeats of MEAN pairwise osc phase difference.
            %
            thing = squeeze(dmm(:,:,t));
            din = diag(thing);
            dinnMm(t,:) = [ mean(din), std(din) ]; 
            % *** dinnMm is MEAN & STD across [num_clust] single clusters of MEAN across simulation 
            %     repeats of MEAN pairwise relative distance between oscillators in SAME cluster
            %
            dout = triu(thing,1);
            dout = dout(dout>0);
            if isempty(dout)
                dout = nan;
            end
            douttMm(t,:) = [ mean(dout), std(dout) ];
            % *** douttMm is MEAN & STD across [num_clust*(num_clust-1)/2] cluster pairs of MEAN across simulation 
            %     repeats of MEAN pairwise relative distance between oscillators in DIFFERENT clusters.\
            %
            

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            % % (2).  MEAN across clusters (in, out or tot network) of STD across simulation repeats of MEAN of pairwise osc phase difference.
            %
            thing = squeeze(dsm(:,:,t));
            din = diag(thing);
            dinnSm(t,:) = [ mean(din) ];
            % *** dinnSm is MEAN across [num_clust] single clusters of STD across simulation 
            %     repeats of MEAN pairwise relative distance between oscillators in SAME cluster.
            %
            dout = triu(thing,1);
            dout = dout(dout>0);
            if isempty(dout)
                dout = nan;
            end
            douttSm(t,:) = [ mean(dout) ];
            % *** douttSm is MEAN across [num_clust*(num_clust-1)/2] cluster pairs of STD across simulation 
            %     repeats of MEAN pairwise relative distance between oscillators in DIFFERENT clusters.
            %
            
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            % % (3).  MEAN across clusters (in, out or tot network) of MEAN across simulation repeats of STD of pairwise osc phase difference.
            %
            thing = squeeze(dms(:,:,t));
            din = diag(thing);
            dinnMs(t,:) = [ mean(din) ];
            % *** dinnMs is MEAN across [num_clust] single clusters of MEAN across simulation 
            %     repeats of STD of pairwise relative distance between oscillators in SAME cluster.
            %
            dout = triu(thing,1);
            dout = dout(dout>0);
            if isempty(dout)
                dout = nan;
            end
            douttMs(t,:) = [ mean(dout) ];
            % *** douttMs is MEAN across [num_clust*(num_clust-1)/2] cluster pairs of of MEAN across simulation 
            %     repeats of STD of pairwise relative distance between oscillators in DIFFERENT clusters.
            %
            
            
        end
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %
        % % (4).  Calculate 3 things above for all oscillator pairs in whole network.
        %
        dtottMm = squeeze(dmm(5,1,:)); 
        % *** dtottMm is MEAN across simulation repeats of MEAN pairwise osc phase difference between ALL OSC PAIRS IN NETWORK
        %
        dtottSm = squeeze(dsm(5,1,:));
        % *** dtottSm is STD across simulation repeats of MEAN pairwise osc phase difference between ALL OSC PAIRS IN NETWORK
        %
        dtottMs = squeeze(dms(5,1,:));
        % *** dtottMs is MEAN across simulation repeats of STD pairwise osc phase difference between ALL OSC PAIRS IN NETWORK
        %
        
        
        % Pack all these stats into d_rpt structure containing each P# and V#
        eval(['d_rpt.','P',num2str(p),'.Mm.din = dinnMm;'])
        eval(['d_rpt.','P',num2str(p),'.Sm.din = dinnSm;'])
        eval(['d_rpt.','P',num2str(p),'.Ms.din = dinnMs;'])
        %
        eval(['d_rpt.','P',num2str(p),'.Mm.dout = douttMm;'])
        eval(['d_rpt.','P',num2str(p),'.Sm.dout = douttSm;'])
        eval(['d_rpt.','P',num2str(p),'.Ms.dout = douttMs;'])
        %
        eval(['d_rpt.','P',num2str(p),'.Mm.dtot = dtottMm;'])
        eval(['d_rpt.','P',num2str(p),'.Sm.dtot = dtottSm;'])
        eval(['d_rpt.','P',num2str(p),'.Ms.dtot = dtottMs;'])
        
        
    
    end
    
    % Add explanations about stats into d_rpt structure containing each P# and V#
    d_rpt.stats = d_sum.stats;
    d_rpt.stats.MnS_types{1} = '(1). Mm = MEAN & STD across clusters (in, out or tot network) of MEAN across simulation repeats of MEAN pairwise osc phase difference.';
    d_rpt.stats.MnS_types{2} = '(2). Sm = MEAN across clusters (in, out or tot network) of STD across simulation repeats of MEAN of pairwise osc phase difference.';
    d_rpt.stats.MnS_types{3} = '(3). Ms = MEAN across clusters (in, out or tot network) of MEAN across simulation repeats of STD of pairwise osc phase difference.';
    %
    d_rpt.stats.MnS_explain{1} = '(1). Mm = MEAN (& STD for IN / OUT) across clusters.';
    d_rpt.stats.MnS_explain{2} = '(2). Sm = STD across simulation repeats.';
    d_rpt.stats.MnS_explain{3} = '(3). Ms = STD of oscillators in within each cluster / cluster pair.';
    
    clear dmm dsm dms din dout thing d_sum                                     % these are all temporary variables.
    clear dinnMm dinnSm dinnMs douttMm douttSm douttMs dtottMm dtottSm dtottMs
    
    
end






plot_pairwise_phase_dists = 1;
%
%
% Plot stats on pairwise oscillator distance (across clusters, simulation repeats & oscillator 
% pairs within clusters) for  different time chunking methods for phase trajectories
if(plot_pairwise_phase_dists)
    xtk = round(linspace(1,num_tg,3));
    for d=1:numel(xtk)
        dtime{d} = num2str(d_rpt.stats.traj_times1(xtk(d),2).*kurParams.tau,2);
    end
%     figure

    
    
    % (1). Plot avg pairwise phase distance within cluster, across clusters
    %      and across entire network for different time chunkings and 
    %      different time points in the simulation.
    %
    for p = 1:num_tchunks
        eval(['thing = d_rpt.P',num2str(p),'.Mm;']) % mean / std cluster, mean sim repeats, mean osc. pairs
        dmeanIa = thing.din(:,1);
        dstdIa = thing.din(:,2); 
        dmeanOa = thing.dout(:,1);
        dstdOa = thing.dout(:,2);
        dmeanTa = thing.dtot;
        dmeanIOsub = thing.dout(:,1) - thing.din(:,1);
        dmeanITsub = thing.dtot(:,1) - thing.din(:,1);
        dmeanIOrat = thing.dout(:,1) ./ thing.din(:,1);
        dmeanITrat = thing.dtot(:,1) ./ thing.din(:,1);
        %
        eval(['thing = d_rpt.P',num2str(p),'.Sm;']) % mean cluster, std sim repeats, mean osc. pairs
        dmeanIb = thing.din;
        dmeanOb = thing.dout;
        dmeanTb = thing.dtot;
        %
        eval(['thing = d_rpt.P',num2str(p),'.Ms;']) % mean cluster, mean sim repeats, std osc. pairs
        dmeanIc = thing.din; 
        dmeanOc = thing.dout;
        dmeanTc = thing.dtot;
        %
        figure
        %
        % % Plot d_in, d_out, d_tot with various error bars.
        %
        subplot(2,2,[1,3]), 
        hold on,
        %
        shadedErrorBar([1:num_tg], dmeanIa, dmeanIc,'r-',0.95)    % STD across osc. pairs in clusters
        shadedErrorBar([1:num_tg], dmeanOa, dmeanOc,'b-',0.95)    % 
        shadedErrorBar([1:num_tg], dmeanTa, dmeanTc,'k-',0.95)    % 
        %
        shadedErrorBar([1:num_tg], dmeanIa, dmeanIb,'g--',0.05)   % STD across simulation repeats
        shadedErrorBar([1:num_tg], dmeanOa, dmeanOb,'c--',0.05)   % 
        shadedErrorBar([1:num_tg], dmeanTa, dmeanTb,'y--',0.05)   % 
        %
        errorbar(1:num_tg, dmeanIa, dstdIa,'ro-','LineWidth',2)   % mean across sim. repeats of STD
        errorbar(1:num_tg, dmeanOa, dstdOa,'bo-','LineWidth',2)   % across clusters & cluster pairs
        plot(1:num_tg, dmeanTa,'ko-','LineWidth',2)  
        grid on
        xlim([ 1 num_tg ])
        set(gca,'Xtick',xtk,'XtickLabel',dtime,'FontSize',16,'FontWeight','Bold')
        xlabel('time point (t)','FontSize',18,'FontWeight','Bold')
        ylabel({['\Theta d_{[\color{red}In\color{black},\color{blue}Out\color{black},\color{black}Tot\color{black}]} / |t|'],...
                ['(1. \sigma across ',num2str(num_clust),' \color{red}clu\color{blue}ste\color{black}r (pair) s)'],...
                ['(2. \sigma across ',num2str(R),' sim. \color{green}re\color{cyan}pe\color{yellow}ats\color{black})'],...
                ['(3. \sigma across osc pairs in cluster [bars])']},...
                'FontSize',18,'FontWeight','Bold') % {'distance per time point','d / |t|'}
        title({'Mean Pairwise \Theta Dist','across clusters & sim. repeats'})
        %  
        % % Plot subtractive difference for (out - in) and (tot - in)
        %
        subplot(2,2,2), hold on, 
        plot(1:num_tg, dmeanIOsub,'b-','LineWidth',2)
        plot(1:num_tg, dmeanITsub,'k-','LineWidth',2)
        xlim([ 1 num_tg ])
        grid on
        ylabel('\Delta\Theta d','FontSize',18,'FontWeight','Bold') % {'distance per time point','d / |t|'}
        set(gca,'Xtick',xtk,'XtickLabel',dtime,'FontSize',16,'FontWeight','Bold')
        title('Difference : d_{[\color{blue}Out-In\color{black}]} & d_{[\color{black}Tot-In\color{black}]}')
        %
        % % Plot subtractive difference for (out - in) and (tot - in)
        %
        subplot(2,2,4), hold on, 
        plot(1:num_tg, dmeanIOrat,'b-','LineWidth',2)
        plot(1:num_tg, dmeanITrat,'k-','LineWidth',2)  
        xlim([ 1 num_tg ])
        grid on
        set(gca,'Xtick',xtk,'XtickLabel',dtime,'FontSize',16,'FontWeight','Bold')
        ylabel('\Delta\Theta d','FontSize',18,'FontWeight','Bold') % {'distance per time point','d / |t|'}
        title('Ratio : d_{[\color{blue}Out/In\color{black}]} & d_{[\color{black}Tot/In\color{black}]}')
        %
        % % Figure title and save
        %
        annotation('textbox', [0 0.9 1 0.1],'String',['Time Chunking ',d_rpt.stats.traj_tpts{p}], ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        %                         hDynDir = ['./output/explore_modularity_via_handcrafted.networks/Figures/Kur_Oscillator_Phase_Dynamics/']
        %                         if(~exist(hDynDir,'dir'))
        %                             mkdir(hDynDir)
        %                         end
        %                         saveGoodImg(hDyn,[hDynDir,methods{m},'_t',num2str(j)],[0 0 1 1])
            

    end
    
end
    
    disp('NOTE: Could also plot average relative pairwise distance stats between oscillators in Velocity.')
    
    
