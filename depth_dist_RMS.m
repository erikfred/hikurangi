% depth_dist_RMS.m
%
% Generate tables that catalog how effective station pair differences are
%

clear; close all

yr_str={'I','III','IV','Vb','VI','VII','VIII','I'};

for i=7:length(yr_str)
    if i<=7
        load(['../processed_data/HOBITSS_' yr_str{i} '_dedrifted_cork.mat'])
    else
        load(['../processed_data/GONDOR_' yr_str{i} '_dedrifted.mat'])
    end

    % sort by depth (makes table more intuitive)
    [stadepth,id]=sort(stadepth);
    stalat=stalat(id);
    stalon=stalon(id);
    staname=staname(id);
    scor=scor(id);
    tf=tf(id);

    count2=0;
    for j=1:length(scor)
        count1=count2+1;
        % station-by-station info extraction
        for k=1:length(scor)
            count2=count2+1;
            
            TBL3{count2,1}=staname{k};
            if isempty(scor{j}) || isempty(scor{k})
                continue
            end
            TBL3{count2,2}=stadepth(k);
            TBL3{count2,3}=stadepth(j)-stadepth(k);
            TBL3{count2,4}=lldistkm([stalat(j) stalon(j)],[stalat(k) stalon(k)]);
            TBL3{count2,5}=rms(scor{k});
            % find intersect of timeseries (accounts for shorter records)
            [~,ia,ib]=intersect(tf{j},tf{k});
            Gp=scor{k}(ib);
            mp=inv(Gp'*Gp)*Gp'*scor{j}(ia);
            TBL3{count2,6}=rms(scor{j}(ia)-mp*scor{k}(ib));
            TBL3{count2,7}=mp;
        end
    end

    % save table
    tbl3=cell2table(TBL3,'VariableNames',{'Name' 'Depth' 'Depth Dif' 'Range' ...
        'RMS' 'Dif RMS' 'Scale Factor'});
    if i<=7
        writetable(tbl3,['../tables/differences/HOBITSS-' yr_str{i} '_sin_diffs.csv'])
    else
        writetable(tbl3,['../tables/differences/GONDOR-' yr_str{i} '_sin_diffs.csv'])
    end

    % plots
    L=length(scor);
    figure(78); clf; hold on
    cmo=cmocean('balance');
    hs=[];
    for jj=1:L
        orig_color=[TBL3{[1:jj-1,jj+1:L]+L*(jj-1),6}]';
        sym_color=floor([TBL3{[1:jj-1,jj+1:L]+L*(jj-1),6}]'*10)/10; % rounds down to 0.1 for discrete color scale
        hs(jj)=scatter([TBL3{[1:jj-1,jj+1:L]+L*(jj-1),3}],...
            [TBL3{[1:jj-1,jj+1:L]+L*(jj-1),4}],250,orig_color,'o','filled');
    end
    % xline(180,'k--')
    % legend(hs,cat(1,leg{:}))
    box on; grid on
    hcb=colorbar;
    xlim([0 2500])
    clim([0.5 1.5])
    colormap(cmo(40:216,:))
    set(gca,'fontsize',16)
    set(hcb,'fontsize',14)
    ylabel(hcb,'RMS (cm)')
    ylabel('Separation (km)')
    xlabel('\DeltaDepth (m)')
    % title('Depth > 300 m')
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 5.5];

    print(['../figures/exploratory/HOBITSS_' yr_str{i} '/differences/depth-dist-rms/sin'],'-dpng','-r300')
    print(['../figures/exploratory/HOBITSS_' yr_str{i} '/differences/depth-dist-rms/sin'],'-depsc','-painters')

    % clear variables
    clearvars('TBL3','tbl3')
end