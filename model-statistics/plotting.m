% Plotting script for GRL, Lu et al. (2026)

% Figure 1-1
% load file
filename  = './model-statistics/statistics-8-13.txt';
formatSpec = '%d %f %f %f %f';% note this assumes 4 columns the last of which is text
fh =fopen(filename,'r');
dataArray = textscan(fh, formatSpec,'HeaderLines',18);
fclose(fh);

fig1 = figure();
fig1.Position(1:4)=[200 150 1000 200];
% find when melt start to produce
ind_start = find((dataArray{3}+dataArray{4})>0,1);
start_time = dataArray{2}(ind_start)/3600/24/365/1e6;
mask = dataArray{2}>0;
mask(1:ind_start) = false;
time = dataArray{2}(mask);
time = time - time(1);

% melt rate (km3/yr)
melt_rate_total = (dataArray{3}(mask)+dataArray{4}(mask))/1e9*3600*24*365;
melt_rate_ecl = dataArray{4}(mask)/1e9*3600*24*365;
% smoothing the curve
melt_rate_total = medfilt1(melt_rate_total, 5);
melt_rate_ecl = medfilt1(melt_rate_ecl, 5);
melt_rate_pyr = melt_rate_total-melt_rate_ecl;
cumulate_melt = zeros(length(time),1);
for j = 2:length(cumulate_melt)
    cumulate_melt(j) = cumulate_melt(j-1)+melt_rate_total(j)*(time(j)-time(j-1))/(3600*24*365);
end
% carbon rate (mol/yr)
cumulate_carbon = zeros(length(time),1);
carbon = dataArray{5}(mask)/1e3*3300/44.01;
for j = 2:length(cumulate_carbon)
    cumulate_carbon(j) = cumulate_carbon(j-1) + carbon(j);
end
% smoothing the curve
carbon_smooth = downsample(cumulate_carbon,10);
time_smooth = downsample(time,10);
carbon_rate = zeros(length(time_smooth),1);
for k = 2:length(carbon_rate)
    carbon_rate(k) = (carbon_smooth(k)-carbon_smooth(k-1))/( (time_smooth(k)-time_smooth(k-1))/365/24/3600 );
end

% plotting
ax1 = gca;
ax1.Position = [0.08,0.18,0.85,0.73];
plot(ax1,time/3600/24/356/1e6+start_time,melt_rate_pyr,'Color','#77AC30','LineWidth',2)
hold on
plot(ax1,time/3600/24/356/1e6+start_time,melt_rate_ecl,'Color','k','LineWidth',1)
hold on
ylim(ax1,[1e-1 100])
xlim(ax1,[112 225])
ax1_pos = ax1.Position;
ax1.YScale = 'log';
xticks(ax1,[112,125,150,175,200,225])
ax1.FontSize = 11;

ax2 = axes('Position',ax1_pos);
plot(ax2,time_smooth/3600/24/356/1e6+start_time,carbon_rate,'Color','#0072BD','LineWidth',1)
ax2.YAxisLocation ='right';
ax2.YColor = '#0072BD';
ax2.Color = 'none';
ax2.YAxis.Scale = 'log';
ax2.XTickLabel = '';
xticks(ax2,[])
ax2.FontSize = 11;
ylim(ax2,[1e11 1e12])
xlim(ax2,[112 225])
ax1.Position = ax2.Position;

ylabel(ax1,'melt rate (km^3/yr)','FontSize',14)
ylabel(ax2,'carbon rate (mol/yr)','FontSize',14)

xlabel(ax1,'time (Ma)','FontSize',14)
xline(ax1,167,'-k','167 Ma','LineWidth',1.5,'Color',"#A2142F",'LabelHorizontalAlignment','right','LabelOrientation','horizontal','FontSize',12)
xline(ax1,187,'-k','187 Ma','LineWidth',1.5,'LabelHorizontalAlignment','right','Color',"#A2142F",'LabelOrientation','horizontal','FontSize',12)
legend(ax1,"peridotite","eclogite",'FontSize',14,'Location','northwest')

%% Change to plot the desired files
filegroup = {
            % 120 km litho + 0-25% Ecl + 350-650 K Tex
            'model-statistics/statistics-8-1.txt';
            'model-statistics/statistics-8-2.txt';
            'model-statistics/statistics-8-3.txt';
            'model-statistics/statistics-8-4.txt';
            'model-statistics/statistics-8-5.txt';
            'model-statistics/statistics-8-6.txt';
            'model-statistics/statistics-8-7.txt';
            'model-statistics/statistics-8-8.txt';
            'model-statistics/statistics-8-9.txt';
            'model-statistics/statistics-8-10.txt';
            'model-statistics/statistics-8-11.txt';
            'model-statistics/statistics-8-12.txt';
            'model-statistics/statistics-8-13.txt';
            'model-statistics/statistics-8-14.txt';
            'model-statistics/statistics-8-15.txt';
            'model-statistics/statistics-8-16.txt';
            'model-statistics/statistics-8-17.txt';
            'model-statistics/statistics-8-18.txt';
            'model-statistics/statistics-8-19.txt';
            'model-statistics/statistics-8-20.txt';

%              % 90 km litho + 0-25% Ecl + 350-650 K Tex
%              'model-statistics/statistics-7-1.txt';
%              'model-statistics/statistics-7-2.txt';
%              'model-statistics/statistics-7-3.txt';
%              'model-statistics/statistics-7-4.txt';
%              'model-statistics/statistics-7-5.txt';
%              'model-statistics/statistics-7-6.txt';
%              'model-statistics/statistics-7-7.txt';
%              'model-statistics/statistics-7-8.txt';
%              'model-statistics/statistics-7-9.txt';
%              'model-statistics/statistics-7-10.txt';
%              'model-statistics/statistics-7-11.txt';
%              'model-statistics/statistics-7-12.txt';
%              'model-statistics/statistics-7-13.txt';
%              'model-statistics/statistics-7-14.txt';
%              'model-statistics/statistics-7-15.txt';             
%              'model-statistics/statistics-7-16.txt';
%              'model-statistics/statistics-7-17.txt';
%              'model-statistics/statistics-7-18.txt';
%              'model-statistics/statistics-7-19.txt';
%              'model-statistics/statistics-7-20.txt';

%              % 60-180 km litho + 15% Ecl + 450 K Tex
%              'model-statistics/statistics-7-33.txt';
%              'model-statistics/statistics-7-13.txt';
%              'model-statistics/statistics-8-13.txt';
%              'model-statistics/statistics-7-32.txt';
%              'model-statistics/statistics-7-31.txt';

             };

formatSpec = '%d %f %f %f %f';% note this assumes 4 columns the last of which is tex
%% Figure 3 and Figure S4
% Compute melt rate and carbon release rate
total_melt5 = zeros(length(filegroup),1);
total_melt10 = zeros(length(filegroup),1);
total_carbon5 = zeros(length(filegroup),1);
total_carbon10 = zeros(length(filegroup),1);
starting_time = zeros(length(filegroup),1);

for i =1:length(filegroup)
    filename  = filegroup{i};
    fh =fopen(filename,'r');
    dataArray = textscan(fh, formatSpec,'HeaderLines',18);
    fclose(fh);
    
    % find when melt start to produce
    ind_start = find((dataArray{3}+dataArray{4})>1,1);
    mask = dataArray{2}>0;
    mask(1:ind_start) = false;
    time = dataArray{2}(mask);
    starting_time(i) = time(1)/3600/24/365/1e6;
    time = time - time(1);
    % melt rate
    melt_rate_total = (dataArray{3}(mask)+dataArray{4}(mask))/1e9*3600*24*365;
    melt_rate_total = medfilt1(melt_rate_total, 5);
    melt_rate_ecl = dataArray{4}(mask)/1e9*3600*24*365;
    melt_rate_ecl = medfilt1(melt_rate_ecl, 5);
    melt_rate_pyr = melt_rate_total-melt_rate_ecl;

    cumulate_melt = zeros(length(time),1);
    for j = 2:length(cumulate_melt)
        cumulate_melt(j) = cumulate_melt(j-1)+melt_rate_total(j)*(time(j)-time(j-1))/(3600*24*365);
    end
    total_melt5(i) = cumulate_melt(find( time>=5e6*3600*24*365 , 1 ));
    total_melt10(i) = cumulate_melt(find( time>=10e6*3600*24*365 , 1 ));

    carbon_rate = zeros(length(time),1);
    cumulate_carbon = zeros(length(time),1);
    carbon = dataArray{5}(mask)/1e3*3300/44.01;
    for j = 2:length(carbon_rate)
        carbon_rate(j) = carbon(j)/( (time(j)-time(j-1))/365/24/3600 );
        cumulate_carbon(j) = cumulate_carbon(j-1) + carbon(j);
    end
    carbon_rate = medfilt1(carbon_rate, 5);
    total_carbon5(i) = cumulate_carbon(find( time>=5e6*3600*24*365 , 1 ));
    total_carbon10(i) = cumulate_carbon(find( time>=10e6*3600*24*365 , 1 ));
end
%% Figure 4-1 and Figure S5
ecl = [5 10 15 20 25];
Tex = [350 450 550 650];
[X,Y] = meshgrid(Tex,ecl);
tempp = linspace(350,650,41);
eclogite = linspace(5,25,21);
[Xq,Yq] = meshgrid(tempp,eclogite);

% scatter with contours
fig1 = figure();
fig1.Position(1:4)=[200 150 620 620];
tiledlayout(2,2,"TileSpacing", "tight");
Z1_raw = reshape(flip(total_melt5),5,4);
Z1 = interp2(X,Y,Z1_raw,Xq,Yq,'cubic');
Z2_raw = reshape(flip(total_melt10),5,4);
Z2 = interp2(X,Y,Z2_raw,Xq,Yq,'cubic');
Z3_raw = reshape(flip(total_carbon5),5,4);
Z3 = interp2(X,Y,Z3_raw,Xq,Yq,'cubic');
Z4_raw = reshape(flip(total_carbon10),5,4);
Z4 = interp2(X,Y,Z4_raw,Xq,Yq,'cubic');
total_num = 20;

nexttile
contour(Xq,Yq,Z1/1e6,[1 10 50 100 200],'k','ShowText','on')
hold on
scatter(reshape(X,total_num,1),reshape(Y,total_num,1),100,flip(total_melt5)/1e6,'filled','MarkerEdgeColor','k')
ylabel("Ecl %",'FontSize',14)
title("Melt Volume 0-5 Myrs",'FontSize',12)
c1 =colorbar();
c1.Ticks = [1 10 50 100 150 200];c1.FontSize = 12;
c1.Limits = [0 200];
clim([0,200])
c1.Title.String = 'x10^6 km^3';c1.Label.FontSize = 12;
colormap(turbo)

nexttile
contour(Xq,Yq,Z2/1e6,[1 10 50 100 200],'k','ShowText','on')
hold on
scatter(reshape(X,total_num,1),reshape(Y,total_num,1),100,flip(total_melt10)/1e6,'filled','MarkerEdgeColor','k')
title("Melt Volume 0-10 Myrs",'FontSize',12)
c2 =colorbar();
c2.Ticks = [1 10 50 100 150 200];c2.FontSize = 12;
clim([0,200])
c2.Limits = [0 200];
c2.Title.String = 'x10^6 km^3';c2.Label.FontSize = 12;

nexttile
contour(Xq,Yq,Z3/1e18,[1 2 3 4 5],'k','ShowText','on')
hold on 
scatter(reshape(X,total_num,1),reshape(Y,total_num,1),100,flip(total_carbon5)/1e18,'filled','MarkerEdgeColor','k')
xlabel("Excess T (K)",'FontSize',14)
title("CO_2 amount 0-5 Myrs",'FontSize',12)
ylabel("Ecl %",'FontSize',14)
c3 =colorbar();
c3.Ticks = [1 3 5 10];c3.FontSize = 12;
clim([0,10])
c3.Limits = [0 10];
c3.Title.String = 'x10^{18} mol';c3.Label.FontSize = 12;

nexttile
contour(Xq,Yq,Z4/1e18,[1 3 5 7 10 15],'k','ShowText','on')
hold on
scatter(reshape(X,total_num,1),reshape(Y,total_num,1),100,flip(total_carbon10)/1e18,'filled','MarkerEdgeColor','k')
title("CO_2 amount 0-10 Myrs",'FontSize',12)
xlabel("Excess T (K)",'FontSize',14)
c4 =colorbar();
c4.Ticks = [1 3 5 10 15 20];c4.FontSize = 12;
clim([0,20])
c4.Limits = [0 20];
c4.Title.String = 'x10^{18} mol';c4.Label.FontSize = 12;

%% Figure 3 and Figure S4
% load average excess temperature from models
aveT_pyr = readmatrix("./model-statistics/aveT_pyr_90.csv");
aveT_ecl = readmatrix("./model-statistics/aveT_ecl_90.csv");
aveP_pyr = readmatrix("./model-statistics/aveP_pyr_90.csv");
aveP_ecl = readmatrix("./model-statistics/aveP_ecl_90.csv");
time_record = readmatrix("./model-statistics/time_90.csv");

aveT_pyr(aveT_pyr == 0) = NaN;
aveT_ecl(aveT_pyr == 0) = NaN;
aveP_pyr(aveT_pyr == 0) = NaN;
aveP_ecl(aveT_pyr == 0) = NaN;

fig1 = figure();
fig1.Position(1:4)=[200 150 1200 900];

for i =1:length(filegroup)
    filename  = filegroup{i};
    fh =fopen(filename,'r');
    dataArray = textscan(fh, formatSpec,'HeaderLines',18);
    fclose(fh);
    
    % find when melt start to produce
    ind_start = find((dataArray{3}+dataArray{4})>0,1);
    mask = dataArray{2}>0;
    mask(1:ind_start) = false;
    start_time = dataArray{2}(ind_start)/3600/24/365/1e6;
    time_T = time_record(i,:) - start_time;
    time_T(time_T < 0) = NaN;
    time = dataArray{2}(mask);
    time = time - time(1);
    melt_rate_total = (dataArray{3}(mask)+dataArray{4}(mask))/1e9*3600*24*365;
    melt_rate_total = medfilt1(melt_rate_total, 5);
    melt_rate_ecl = dataArray{4}(mask)/1e9*3600*24*365;
    melt_rate_ecl = medfilt1(melt_rate_ecl, 5);
    melt_rate_pyr = melt_rate_total-melt_rate_ecl;

    carbon_rate = zeros(length(time),1);
    cumulate_carbon = zeros(length(time),1);
    carbon = dataArray{5}(mask)/1e3*3300/44.01;
    for j = 2:length(carbon_rate)
        cumulate_carbon(j) = cumulate_carbon(j-1) + carbon(j);
    end
    for k = 2:length(carbon_rate)
        carbon_rate(k) = (cumulate_carbon(k)-cumulate_carbon(k-1))/( (time(k)-time(k-1))/365/24/3600 );
    end
    carbon_rate = medfilt1(carbon_rate, 10);

    % plot a-d
    subplot(4,3,i*3-2);
    plot(time/3600/24/356/1e6,hhh,'LineWidth',1,'Color',"#0072BD")
    hold on
    plot(time/3600/24/356/1e6,melt_rate_ecl,'LineWidth',1,'Color',"#D95319")
    ax1 = gca;
    ax1.YAxis.Scale = 'log';
    hor_loc = 1-i*0.235;
    ax1.Position = [0.04 hor_loc 0.267 0.205];
    yticks(ax1,[1e-3 1e-1 1 10])
    ylabel(ax1,'melt rate (km^3/yr)','FontSize',14)
    xlabel(ax1,'time (Ma)','FontSize',14)
    ylim(ax1,[1e-1 50])
    xlim(ax1,[0 60])
    if i == 4
        legend(ax1,"peridotite","eclogite",Orientation="horizontal",Location="southeast",fontsize=14)
        ylim(ax1,[5e-3 10])
    end
    ax1.FontSize = 12;
    xlabel(ax1,'time (Ma)','FontSize',14)

    % plot e-h
    subplot(4,3,i*3-1);
    ax2 = gca;
    ax2.Position = [0.36 hor_loc 0.269 0.205];
    ax2_pos = ax2.Position;
    colororder({'k','k'})
    yyaxis left
    plot(ax2,time_T,aveT_pyr(i,:),'Marker','o','LineStyle','-','LineWidth',1,'Color',"#0072BD")
    hold on
    plot(ax2,time_T,aveT_ecl(i,:),'Marker','square','LineStyle','-','LineWidth',1,'Color',"#D95319")
    hold on
    ylim(ax2,[100 400])
    if i == 4
        ylim(ax2,[50 350])
    end
    yticks(ax2,[100 200 300 400])
    ylabel(ax2,'T_{Excess} (K)','FontSize',14)

    yyaxis right
    plot(ax2,time_T,aveP_pyr(i,:)/1e9,'Marker','o','LineStyle',':','LineWidth',1,'Color',"#0072BD")
    hold on
    plot(ax2,time_T,aveP_ecl(i,:)/1e9,'Marker','square','LineStyle',':','LineWidth',1,'Color',"#D95319")
    ylim(ax2,[2 12])
    yticks(ax2,[4 6 8 10])
    ylabel(ax2,'Pressure (GPa)','FontSize',14)
    ax2.FontSize = 12;
    xlabel(ax2,'time (Ma)','FontSize',14)
    if i == 4
       ylim(ax2,[0 10])
       legend(ax2,"peri-\DeltaT","ecl-\DeltaT","peri-P","ecl-P",Orientation="horizontal",Location="northeast",fontsize=14)
    end    
    xlim(ax2,[0 60])

    % plot i-l
    subplot(4,3,i*3);
    plot(time/3600/24/356/1e6,carbon_rate,'LineWidth',1,'Color',"#D95319")
    ax3 = gca;
    ax3.YAxis.Scale = 'log';
    ax3.Position = [0.725 hor_loc 0.267 0.205];
    yticks(ax3,[1e10 1e11 1e12])
    ylabel(ax3,'CO_2 rate (mol/yr)','FontSize',14)
    xlabel(ax3,'time (Ma)','FontSize',14)
    if i == 4
        ylim(ax3,[2e10 1e12])
        yline(ax3,0.8e11,'--k','Tail phase value','LineWidth',1.5,'LabelHorizontalAlignment','left');
    else
        ylim(ax3,[8e10 3e12])
        yline(ax3,1.2e11,'--k','Tail phase value','LineWidth',1.5,'LabelHorizontalAlignment','left');
    end
    xlim(ax3,[0 60])
    ax3.FontSize = 12;
    xlabel(ax3,'time (Ma)','FontSize',14)    
end

%% Figure S3
fig1 = figure();
fig1.Position(1:4)=[200 150 400 900];
color_list = {[0.95 0 0];[0.850 0.2 0.0080];[0.8 0.3 0];[0.8000 0.5350 0.0080];[0.9290 0.6940 0.3250]};
for i =1:length(filegroup)
    filename  = filegroup{i};
    fh =fopen(filename,'r');
    dataArray = textscan(fh, formatSpec,'HeaderLines',18);
    fclose(fh);
    
    % find when melt start to produce
    ind_start = find((dataArray{3}+dataArray{4})>1,1);
    mask = dataArray{2}>0;
    mask(1:ind_start) = false;
    time = dataArray{2}(mask);
    starting_time(i) = time(1)/3600/24/365/1e6;
    time = time - time(1);
    % melt rate
    melt_rate_total = (dataArray{3}(mask)+dataArray{4}(mask))/1e9*3600*24*365;
    melt_rate_total = medfilt1(melt_rate_total, 5);
    melt_rate_ecl = dataArray{4}(mask)/1e9*3600*24*365;
    melt_rate_ecl = medfilt1(melt_rate_ecl, 5);
    melt_rate_pyr = melt_rate_total-melt_rate_ecl;

    cumulate_melt = zeros(length(time),1);
    for j = 2:length(cumulate_melt)
        cumulate_melt(j) = cumulate_melt(j-1)+melt_rate_total(j)*(time(j)-time(j-1))/(3600*24*365);
    end

    carbon_rate = zeros(length(time),1);
    cumulate_carbon = zeros(length(time),1);
    carbon = dataArray{5}(mask)/1e3*3300/44.01;
    % carbon rate
    for j = 2:length(carbon_rate)
        cumulate_carbon(j) = cumulate_carbon(j-1) + carbon(j);
    end
    for k = 2:length(carbon_rate)
        carbon_rate(k) = (cumulate_carbon(k)-cumulate_carbon(k-1))/( (time(k)-time(k-1))/365/24/3600 );
    end
    carbon_rate = medfilt1(carbon_rate, 10);

    subplot(4,1,1);
    plot(time/3600/24/356/1e6,melt_rate_total,'LineWidth',1,'Color',color_list{i})
    hold on
    ax1 = gca;
    
    subplot(4,1,2);
    plot(time/3600/24/356/1e6,cumulate_melt,'LineWidth',2,'Color',color_list{i})
    hold on
    ax2 = gca;

    subplot(4,1,3);
    plot(time/3600/24/356/1e6,carbon_rate,'LineWidth',1,'Color',color_list{i})
    hold on
    ax3 = gca;

    subplot(4,1,4);
    plot(time/3600/24/356/1e6,cumulate_carbon,'LineWidth',2,'Color',color_list{i})
    hold on
    ax4 = gca;
end

ax1.YAxis.Scale = 'log';
ax3.YAxis.Scale = 'log';

ax1.Position = [0.15 0.72 0.8 0.27];
ax2.Position = [0.15 0.52 0.8 0.14];
ax3.Position = [0.15 0.24 0.8 0.25];
ax4.Position = [0.15 0.05 0.8 0.14];

ax1.FontSize = 12;ax2.FontSize = 12;ax3.FontSize = 12;ax4.FontSize = 12;

yticks(ax1,[1e-3 1e-1 1 10])
yticks(ax1,[1e-3 1e-1 1 10])

ylabel(ax1,'Total melt rate (km^3/yr)','FontSize',14)
ylabel(ax2,'cumulative melt (km^3)','FontSize',14)
ylabel(ax3,'carbon rate (mol/yr)','FontSize',14)
ylabel(ax4,'cumulative carbon (mol)','FontSize',14)

xlabel('time (Ma)','FontSize',14)
% legend(ax1,"650K","550K","450K","350K",Orientation="horizontal",Location="southeast",fontsize=14)
% legend(ax2,"180km","150km","120km","90km","60km",Orientation="horizontal",Location="northwest",fontsize=14)
legend(ax2,"25%","20%","15%","10%","5%",Orientation="horizontal",Location="northeast",fontsize=14)


ylim(ax1,[1e-3 10])
ylim(ax2,[0 2e8])
ylim(ax3,[8e10 3e12])
ylim(ax4,[0 2e19])
xlim(ax1,[0 60])
xlim(ax2,[0 60])
xlim(ax3,[0 60])
xlim(ax4,[0 60])
