% difference_all.m
%
% Difference all possible combinations of stations, after various versions
% of drift correction
%

clear; close all

yr_str={'I','III','IV','Vb','VI','VII','VIII'};

for i=1:length(yr_str)
    load(['../processed_data/HOBITSS_' yr_str{i} '_dedrifted.mat'])
    
    if i==3 % wellheads come online
        cork = load('../processed_data/CORK_dedrifted.mat');

        stadepth=cat(2,stadepth,cork.stadepth(2:3));
        stalat=cat(2,stalat,cork.stalat(2:3));
        stalon=cat(2,stalon,cork.stalon(2:3));
        staname=cat(2,staname,cork.staname(2:3));

        % determine temporal overlap
        t_min=min([tf{:}]);
        t_max=max([tf{:}]);
        cond=(cork.tf{2}>=t_min) & (cork.tf{2}<=t_max);
        tf=cat(2,tf,{cork.tf{2}(cond)});
        pcor=cat(2,pcor,{[]});
        scor=cat(2,scor,{cork.scor{2}(cond)});
        lcor=cat(2,lcor,{cork.lcor{2}(cond)});
        cond=(cork.tf{3}>=t_min) & (cork.tf{3}<=t_max);
        tf=cat(2,tf,{cork.tf{3}(cond)});
        pcor=cat(2,pcor,{[]});
        scor=cat(2,scor,{cork.scor{3}(cond)'});
        lcor=cat(2,lcor,{cork.lcor{3}(cond)'});
    elseif i>3 % wellheads + BPR online
        cork = load('../processed_data/CORK_dedrifted.mat');
        stadepth=cat(2,stadepth,cork.stadepth(1:3));
        stalat=cat(2,stalat,cork.stalat(1:3));
        stalon=cat(2,stalon,cork.stalon(1:3));
        staname=cat(2,staname,cork.staname(1:3));
        
        % determine temporal overlap
        t_min=min([tf{:}]);
        t_max=max([tf{:}]);
        cond=(cork.tf{1}>=t_min) & (cork.tf{1}<=t_max);
        tf=cat(2,tf,{cork.tf{1}(cond)});
        pcor=cat(2,pcor,{[]});
        scor=cat(2,scor,{cork.scor{1}(cond)});
        lcor=cat(2,lcor,{cork.lcor{1}(cond)});
        cond=(cork.tf{2}>=t_min) & (cork.tf{2}<=t_max);
        tf=cat(2,tf,{cork.tf{2}(cond)});
        pcor=cat(2,pcor,{[]});
        scor=cat(2,scor,{cork.scor{2}(cond)});
        lcor=cat(2,lcor,{cork.lcor{2}(cond)});
        cond=(cork.tf{3}>=t_min) & (cork.tf{3}<=t_max);
        tf=cat(2,tf,{cork.tf{3}(cond)});
        pcor=cat(2,pcor,{[]});
        scor=cat(2,scor,{cork.scor{3}(cond)});
        lcor=cat(2,lcor,{cork.lcor{3}(cond)});
    end
    
    for l=1:3
        % sort by depth
        [d,id] = sort(stadepth,'descend');
        t = tf(id);
        lat = stalat(id);
        lon = stalon(id);
        sta = staname(id);

        if l==1
            p=pcor(id);
        elseif l==2
            p=scor(id);
        else
            p=lcor(id);
        end

        % remove empty fields, if any
        t=t(~cellfun(@isempty,p));
        d=d(~cellfun(@isempty,p));
        lat=lat(~cellfun(@isempty,p));
        lon=lon(~cellfun(@isempty,p));
        sta=sta(~cellfun(@isempty,p));
        p=p(~cellfun(@isempty,p));

        for k=1:length(p)
            t1 = t{k};
            p1 = p{k};
            if length(t1)>length(p1)
                t1=t1(1:length(p1));
            end
            for j=1:length(p)
                t2 = t{j};
                p2 = p{j};
                if length(t2)>length(p2)
                    t2=t2(1:length(p2));
                end

                [~,ia,ib] = intersect(round(t1,6),round(t2,6));
                tdif{k,j} = t1(ia);
                pdif{k,j} = p2(ib) - p1(ia);
                xydif(j) = round(lldistkm([lat(k),lon(k)],[lat(j),lon(j)]));
            end

            %-----PLOTTING-----%
            if ~exist(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/lin/HOBITSS_' yr_str{i} '_' sta{k}],'dir')
                mkdir(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/poly/HOBITSS_' yr_str{i} '_' sta{k}]);
                mkdir(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/sin/HOBITSS_' yr_str{i} '_' sta{k}]);
                mkdir(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/lin/HOBITSS_' yr_str{i} '_' sta{k}]);
            end
            figure(11); clf; hold on
            figure(12); clf;
            for j=1:length(p)
                if isempty(pdif{k,j})
                    continue
                end
                figure(11) % stack fig
                plot(tdif{k,j},pdif{k,j}+5*(j-1),'linewidth',1)
                text(tdif{k,j}(end)+10,pdif{k,j}(end)+5*(j-1),{[sta{j} ' - ' sta{k}];...
                    ['\Deltad = ' num2str(d(j)-d(k)) ' m, \Deltaxy = '...
                    num2str(xydif(j)) ' km']},'fontsize',14)

                figure(12) % individual fig
                plot(tdif{k,j},pdif{k,j},'linewidth',1)
                text(tdif{k,j}(1),4,['RMS = ' num2str(std(pdif{k,j})) ' cm'],'fontsize',12)
                title([sta{j} ' - ' sta{k} ', \Deltad = ' num2str(d(j)-d(k)) ...
                    ' m, \Deltaxy = ' num2str(xydif(j)) ' km'],'fontsize',14)

                xlim([min(cat(2,tf{:}))-30 max(cat(2,tf{:}))+30])
                datetick('x','keeplimits')
                ylim([-10 10])
                set(gca,'fontsize',14)
                ylabel('\DeltaP (hPa)')
                box on; grid on

                fh=gcf;
                fh.PaperUnits='inches';
                fh.PaperPosition=[0 0 11 8.5];
                if l==1
                    print(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                        '/differences/poly/HOBITSS_' yr_str{i} '_' sta{k} '/'...
                        sta{j} '-' sta{k}],'-dpng','-r300')
                elseif l==2
                    print(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                        '/differences/sin/HOBITSS_' yr_str{i} '_' sta{k} '/'...
                        sta{j} '-' sta{k}],'-dpng','-r300')
                else
                    print(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                        '/differences/lin/HOBITSS_' yr_str{i} '_' sta{k} '/'...
                        sta{j} '-' sta{k}],'-dpng','-r300')
                end
            end
            figure(11)
            datetick('x')
            lim_x = xlim;
            lim_x(2) = lim_x(2)+60;
            xlim([lim_x(1) lim_x(2)])
            datetick('x','keeplimits')
            set(gca,'fontsize',14)
            ylabel('\DeltaP (hPa)')
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 8.5 11];
            if l==1
                print(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/poly/HOBITSS_' yr_str{i} '_difstack_' sta{k}],'-dpng','-r300')
            elseif l==2
                print(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/sin/HOBITSS_' yr_str{i} '_difstack_' sta{k}],'-dpng','-r300')
            else
                print(['../figures/exploratory/HOBITSS_' yr_str{i} ...
                    '/differences/lin/HOBITSS_' yr_str{i} '_difstack_' sta{k}],'-dpng','-r300')
            end
        end
        clearvars('tdif','pdif')
    end
end