function [] = S_AltiFunc_CS2(path, DamName, Lv_range, Area, First_filter, Second_filter, Third_filter, low_con, Fourth_filter, runoff_rate, output, RGpath, shppath, numWorkers)

%format long
%addpath(genpath('/path/to/mappingtoolbox'))

numWorkers = numWorkers;

%existingPool = gcp('nocreate');
%    if ~isempty(existingPool)
%    delete(existingPool);
%    end
%    pool = parpool(numWorkers);

fprintf('\nCryoSat-2 analysis starts \n');
fprintf('First_filter = %s, Second_filter = %s, Third_filter = %s, Fourth_filter = %s\n', First_filter, Second_filter, Third_filter, Fourth_filter);
fprintf('\nAnalysing how many tracks passing over water surface when using maximum water mask \n');

path = path;

files = dir(fullfile(path, 'CS_*.nc'));
filename = {files.name};
file_count = length(files);

Passfilename = cell(0);

%DamNames = {'JA', 'SJ', 'DC', 'NG', 'HC', 'IH', 'AD', 'CJ', 'SY'};
SatName = 'CS2';
DamName = DamName;

Lv_range = Lv_range;
Area = Area;

First_filter = First_filter;
Second_filter = Second_filter; % Select First_filter or 'N' 
%interval = 1; %parfor = 1
Third_filter = Third_filter;
low_con = low_con;
Fourth_filter = Fourth_filter;
runoff_rate = runoff_rate; % if the rate 70%, it should be expressed by runoff_rate

output = output;
%mkdir(['',output,''])

shppath = shppath;
Shp_Stream = shaperead(fullfile(shppath, sprintf(('%s_%s.shp'), DamName, num2str(max(Lv_range)))));

BefLv_X = Shp_Stream.X;
N_BefLv_X = BefLv_X(~isnan(BefLv_X));
BefLv_Y = Shp_Stream.Y;
N_BefLv_Y = BefLv_Y(~isnan(BefLv_X));
%P = [reshape(x_Up, [],1), reshape(y_Up,  [],1)];
%PS_Up = polyshape(P, 'Simplify', false);

mean_x_Up = mean(N_BefLv_X);

utmZone = floor((mean_x_Up + 180)/6) + 1;
if utmZone > 60
    utmZone = 60;
end

utmCrs = projcrs(32600 + utmZone, 'Authority', 'EPSG');

l = 1;

Sel_Set_FL = [];
Passfilename = [];

epoch = datetime(2000,1,1,0,0,0,'TimeZone','UTC');

RGpath = RGpath;
%Real_Gauged = readtable(fullfile(RGpath, ['Real' DamName 'Level.xlsx']));
%Gauged_time = datetime(Real_Gauged.Var1,'TimeZone','UTC');

Real_Gauged_daily = readtable(fullfile(RGpath, ['DailyRain_' DamName '.xlsx']));
Gauged_time_daily = datetime(Real_Gauged_daily.Var1,'TimeZone','UTC');
Gauged_rainfall_daily = Real_Gauged_daily.Var2;
Gauged_daily = table(Gauged_time_daily, Gauged_rainfall_daily);

%Sum_filename = ['!Summary_' output '.csv'];
%csv_path = fullfile(output, Sum_filename);
%headers1 = {'Time','Count','Min','Max','Count_50', 'waterlevel'};


for j = 1 : file_count

% 1_Read an Altimetry NC.file and Make One-Set file
input_file = fullfile(path, filename{j});

lat = ncread(input_file,'lat_20_ku');
lon = ncread(input_file,'lon_20_ku');

% 3_Extract track points passing over the upstream
in = inpolygon(lon, lat, BefLv_X, BefLv_Y);
x_in = lon(in);
y_in = lat(in);

    fprintf(' %d / %d completed', j, file_count);
    pause(0.06);

    if isempty(x_in)
        clear('lon', 'lat');
    else
        time = ncread(input_file,'time_20_ku');
        alt = ncread(input_file,'alt_20_ku')*1.000000000000000e-03;
        range_ocog = ncread(input_file,'range_ocog_20_ku')*1.000000000000000e-03;
        range_ocean = ncread(input_file,'range_ocean_20_ku')*1.000000000000000e-03;
    %iono_cor = ncread(input_file,'iono_cor_alt_20_ku');
    %iono_cor(isnan(iono_cor)) = 0;
    %idx = isnan(WL_CS2);
    %WL_CS2 = WL_CS2(~idx);
    %time = time(~idx);
    %lon = lon(~idx);
    %lat = lat(~idx);
    %W_Type = W_Type(~idx);
    %W_Size = W_Size(~idx);
    %W_source = W_source(~idx);
    %W_ID = W_ID(~idx);
        Set = table(time, lon, lat, alt, range_ocog, range_ocean);
        Find_Set = find(ismember(Set.lon, x_in));
        Sel_Set = Set(Find_Set,:);
    
    %latHZ1=ncread(input_file,'lat_01');
        lonHZ1=ncread(input_file,'lon_01');
    
        neg_idx = lonHZ1 < 0;
        lonHZ1(neg_idx) = [];
        [lonHZ1, idx] = unique(lonHZ1);

        Geoid=ncread(input_file,'geoid_01');
        Geoid(neg_idx) = [];
        Geoid = Geoid(idx);
        sel_Geoid = interp1(lonHZ1, Geoid, Sel_Set.lon);
    
        iono_cor_01 = ncread(input_file,'iono_cor_gim_01');
        iono_cor_01(neg_idx) = [];
        iono_cor_01 = iono_cor_01(idx);
        sel_iono_cor = interp1(lonHZ1, iono_cor_01, Sel_Set.lon);
        sel_iono_cor(isnan(sel_iono_cor)) = 0;

        mod_dry_tropo_cor_01=ncread(input_file,'mod_dry_tropo_cor_01');
        mod_dry_tropo_cor_01(neg_idx) = [];
        mod_dry_tropo_cor_01 = mod_dry_tropo_cor_01(idx);
        sel_dry_tropo_cor = interp1(lonHZ1, mod_dry_tropo_cor_01, Sel_Set.lon);

        mod_wet_tropo_cor_01=ncread(input_file,'mod_wet_tropo_cor_01');
        mod_wet_tropo_cor_01(neg_idx) = [];
        mod_wet_tropo_cor_01 = mod_wet_tropo_cor_01(idx);
        sel_wet_tropo_cor = interp1(lonHZ1, mod_wet_tropo_cor_01, Sel_Set.lon);

        pole_tide_01=ncread(input_file,'pole_tide_01');
        pole_tide_01(neg_idx) = [];
        pole_tide_01 = pole_tide_01(idx);
        sel_pole_tide_cor = interp1(lonHZ1, pole_tide_01, Sel_Set.lon);

        solid_earth_tide_01=ncread(input_file,'solid_earth_tide_01');
        solid_earth_tide_01(neg_idx) = [];
        solid_earth_tide_01 = solid_earth_tide_01(idx);
        sel_earth_tide_cor = interp1(lonHZ1, solid_earth_tide_01, Sel_Set.lon);
    
        Result_ocog = (Sel_Set.alt - Sel_Set.range_ocog)*1000 - sel_Geoid - sel_dry_tropo_cor - sel_wet_tropo_cor - sel_pole_tide_cor - sel_earth_tide_cor - sel_iono_cor;
        Result_ocean = (Sel_Set.alt - Sel_Set.range_ocean)*1000 - sel_Geoid - sel_dry_tropo_cor - sel_wet_tropo_cor - sel_pole_tide_cor - sel_earth_tide_cor - sel_iono_cor;
    
        Sel_Set2 = table(Sel_Set.time, Sel_Set.lon, Sel_Set.lat, Result_ocog);
        Sel_Set3 = table(Sel_Set.time, Sel_Set.lon, Sel_Set.lat, Result_ocean);

        Sel_Set_FL = [Sel_Set_FL; Sel_Set2; Sel_Set3];
    %Sel_Set_FL = [Sel_Set_FL; Sel_Set2];
    
        Passfilename{end+1} = input_file;
        
            if numel(Passfilename) > 1 && strcmp(Passfilename{end}, Passfilename{end-1})
                Passfilename(end) = [];
            end
    end

        fprintf(repmat('\b', 1, numel(sprintf('%d / %d completed', j, file_count)) +1 ));
end

fprintf('The No. of Tracks passing over Max. W.M.: %d of %d \n', size(Passfilename, 2), file_count);

if isempty(Passfilename)
    warning('There is no pass of %s', SatName);
    return;
end

time_KOR = seconds(Sel_Set_FL.Var1) + epoch + hours(9);
Sel_Set_FL = table(time_KOR(:,1), time_KOR(:,1), Sel_Set_FL.Var2, Sel_Set_FL.Var3, Sel_Set_FL.Var4);
Sel_Set_FL.Properties.VariableNames = {'Date', 'Date_sec', 'lon', 'lat', 'WL_CS2'};

save(fullfile(output, [output 'CS2_f0.mat']));


%% Filter 1 : Standard Operation Level (No option - Mandatory)

if strcmpi(First_filter, 'Y')
    fprintf('\nStart Filter 1 - Standard Operation Level \n');
    idx = Sel_Set_FL.WL_CS2 < min(Lv_range) | Sel_Set_FL.WL_CS2 > max(Lv_range);    
    Sel_Set_FL = Sel_Set_FL(~idx, :);
    idx = isnan(Sel_Set_FL.WL_CS2);
    Sel_Set_FL = Sel_Set_FL(~idx, :);
    elseif strcmpi(First_filter, 'N')
        fprintf('\nSkip Filter 1 - Standard Operation Level \n');
        idx = isnan(Sel_Set_FL.WL_CS2);
        Sel_Set_FL = Sel_Set_FL(~idx, :);
end
    
    save(fullfile(output, [output 'CS2_f1.mat']));

%end

%% Filter 2 : Applying Dynamic Water Masks
if strcmpi(Second_filter, 'Y')
    fprintf('Start Filter 2 - Dynamic Water Masks \n');

% Change sec to day for the Column 'Date'
Sel_Set_FL.Date = dateshift(Sel_Set_FL.Date, 'start', 'day');
groupedTable = groupsummary(Sel_Set_FL, 'Date');

clear('x_in', 'x_Up', 'y_in', 'y_Up', 'in');
%fprintf('The No. of Tracks passing over Max. W.M. %d of %d \n', size(Passfilename, 2), file_count);

WaterMasks = [num2cell(Lv_range)];

DyWMs = [];

for e = 1:height(groupedTable)
    current_group_date = groupedTable.Date(e);
    % Select rows corresponding a such date group
    current_rows = Sel_Set_FL.Date == current_group_date;
    Sel_Set_cur_FL = Sel_Set_FL(current_rows, :);
    
    t_test = [];
        
        fprintf('%d / %d completed', e, height(groupedTable))

    parfor w = 1 : numel(WaterMasks)
        WaterMask = WaterMasks{w};
        
Shp_Stream = shaperead(fullfile(shppath, sprintf(('%s_%s.shp'), DamName, num2str(WaterMask))));        

x_Up = Shp_Stream.X;
y_Up = Shp_Stream.Y;

in = inpolygon(Sel_Set_cur_FL.lon, Sel_Set_cur_FL.lat, x_Up, y_Up);
x_in = Sel_Set_cur_FL.lon(in);
y_in = Sel_Set_cur_FL.lat(in);

x_in_table = table(x_in);
find_Set_Dy = find(ismember(Sel_Set_cur_FL.lon, x_in_table.x_in));
Sel_Set_Dy = Sel_Set_cur_FL(find_Set_Dy, :);

Sel_Set_end{w} = Sel_Set_Dy;

[~, p_value_left, ~, stats_left] = ttest2(Sel_Set_Dy.WL_CS2, Sel_Set_cur_FL.WL_CS2, 'tail', 'left');
[~, p_value_right, ~, stats_right] = ttest2(Sel_Set_Dy.WL_CS2, Sel_Set_cur_FL.WL_CS2, 'tail', 'right');
[~, p_value_both, ~, stats_both] = ttest2(Sel_Set_Dy.WL_CS2, Sel_Set_cur_FL.WL_CS2, 'tail', 'both');

t_test_end{w} = [WaterMask, p_value_left, stats_left.tstat, stats_left.df, p_value_right, stats_right.tstat, stats_right.df, p_value_both, stats_both.tstat, stats_both.df];

%clear('x_in', 'x_Up', 'y_in', 'y_Up', 'in');

    end
    
for y = 1 : numel(WaterMasks)
    Sel_Set_Dy = Sel_Set_end{y};
    varName1 = ['Sel_Set_Dy_', datestr(current_group_date, 'yy_mm_dd'), '_', num2str(WaterMasks{y})];
    eval([varName1 ' = Sel_Set_Dy;']);
end

t_test = vertcat(t_test_end{:});
t_test = table(t_test(:,1), t_test(:,2), t_test(:,3), t_test(:,4), t_test(:,5), t_test(:,6), t_test(:,7),  t_test(:,8), t_test(:,9), t_test(:,10));
t_test.Properties.VariableNames = {'WaterMask', 'p_value_left', 't_stat_left', 'df_left', 'p_value_right', 't_stat_right', 'df_right', 'p_value_both', 't_stat_both', 'df_both'};
varName2 = ['t_test', '_',datestr(current_group_date, 'yy_mm_dd')];
%assignin('base', varName2, t_test);
eval([varName2 ' = t_test;']);

minValue = min(abs(eval([varName2 '.p_value_left'])));

if  isnan(minValue)
    DSel_WM = 'FL';
    else
    %mVIndex = evalin('base', [varName2 '.p_value_right == minValue']);
    mVIndex = eval([varName2 '.p_value_left == minValue']);
    %assignin('base', 'mVIndex', mVIndex);
    DSel_WM = min(eval([varName2 '{mVIndex, 1}']));
end

DyWMs = [DyWMs {DSel_WM}];

fprintf(repmat('\b', 1, numel(sprintf('%d / %d completed', e, height(groupedTable))) ));

end

allVariables = who;
Result_al = [];
%T2 = [];

for v = 1 : numel(DyWMs)
        DyWM = DyWMs{v};
        DyWM_date = groupedTable{v, 1};
        varName3 = ['Sel_Set_Dy_', datestr(DyWM_date, 'yy_mm_dd'), '_', num2str(DyWM)];        
        
        if DyWM == 'FL'
            current_rows = Sel_Set_FL.Date == DyWM_date;
            DyWM_Sel_Set = Sel_Set_FL(current_rows, :);
        else
            DyWM_Sel_Set = eval(varName3);
        end

            %Close_Time = find(abs(Gauged_time - datetime(DyWM_Sel_Set{1,2})) < days(1/288));
            %Real_WL = Real_Gauged(Close_Time,:);
            %Real_WL = repmat(Real_WL, size(DyWM_Sel_Set, 1), 1);
            %Result_OHS = DyWM_Sel_Set.WL_CS2 - Real_WL.Var2;
            STDV_track = std(DyWM_Sel_Set.WL_CS2);
            CI = STDV_track / sqrt(length(DyWM_Sel_Set.WL_CS2));
            CI = repmat(CI, size(DyWM_Sel_Set.WL_CS2));

            Result = table(DyWM_Sel_Set{:,1}, DyWM_Sel_Set{:,2}, DyWM_Sel_Set{:,3}, DyWM_Sel_Set{:,4}, DyWM_Sel_Set{:,5}, CI);
            Result.Properties.VariableNames = {'Date', 'Date_sec', 'lon', 'lat', 'WL_CS2', 'CI'};
            
            Result_al = [Result_al; Result];
            %writetable(Result, fullfile(output, sprintf('IS2_%s.xls', varName3)), 'fileType', 'spreadsheet');

            %OHS_min = min(Result_OHS);
            %OHS_max = max(Result_OHS);
            %count_50 = sum(abs(Result_OHS) <= 0.3);
            %waterlevel = Result.gauged_WL(1,1);
            
            %DyWM_date_str = datestr(DyWM_date, 'yy-mm-dd'); % Or any format you like
            %T = table({DyWM_date_str}, size(DyWM_Sel_Set, 1), OHS_min, OHS_max, count_50, waterlevel, 'VariableNames', {'Time', 'Count', 'Min', 'Max', 'Count_50', 'waterlevel'});
            %T2 = [T2; T];
end

%writetable(T2, Sum_filename);

    idx = isnan(Result_al.WL_CS2);
    Result_al = Result_al(~idx, :);

%    WL_Alt_al = Result_al.WL_CS2;
%    WL_Real_al = Result_al.gauged_WL;

save(fullfile(output, [output 'CS2_f2.mat']));
end


%%
if strcmpi(Second_filter, 'N')
    fprintf('Skip Filter 2 - Dynamic Water Masks \n');

    % Change sec to day for the Column 'Date'
    Sel_Set_FL.Date = dateshift(Sel_Set_FL.Date, 'start', 'day');
    groupedTable = groupsummary(Sel_Set_FL, 'Date');

    clear('x_in', 'x_Up', 'y_in', 'y_Up', 'in');
    %fprintf('The No. of Tracks passing over Max. W.M. %d of %d', size(Passfilename, 2), file_count);
    
    DyWM_Sel_Set = Sel_Set_FL;
    
    %for v = 1 : height(DyWM_Sel_Set)
    %Close_Time = find(abs(Gauged_time - datetime(DyWM_Sel_Set{v,2})) < days(1/288));
    %Real_WL(v, :) = Real_Gauged(Close_Time, :);
    %end
    
    %Real_WL = repmat(Real_WL, size(DyWM_Sel_Set, 1), 1);
    %Result_OHS = DyWM_Sel_Set.WL_CS2 - Real_WL.Var2;
    
        CI = zeros(height(DyWM_Sel_Set), 1);
    
    for i = 1:height(groupedTable)
        groupIndices = DyWM_Sel_Set.Date == groupedTable.Date(i);
        STDV_track = std(DyWM_Sel_Set.WL_CS2(groupIndices));
        CI_value = STDV_track / sqrt(sum(groupIndices));
        CI(groupIndices) = CI_value;
    end

    Result_al = table(DyWM_Sel_Set{:,1}, DyWM_Sel_Set{:,2}, DyWM_Sel_Set{:,3}, DyWM_Sel_Set{:,4}, DyWM_Sel_Set{:,5}, CI);
    Result_al.Properties.VariableNames = {'Date', 'Date_sec', 'lon', 'lat', 'WL_CS2', 'CI'};
    
    %Result_al = table(DyWM_Sel_Set{:,1}, DyWM_Sel_Set{:,2}, DyWM_Sel_Set{:,3}, DyWM_Sel_Set{:,4}, DyWM_Sel_Set{:,5}, Result_OHS, Real_WL.Var1, Real_WL.Var2);
    %Result_al.Properties.VariableNames = {'Date', 'Date_sec', 'lon', 'lat', 'WL_CS2', 'WL_CS2_WL', 'Date_gauged', 'gauged_WL'};
    
end

%% Filter 3 : Applying Tracks with lower confidence in data quantity

if strcmpi(Third_filter, 'Y')
    fprintf('Start Filter 3 - Data Reliability \n');
    groupedTable2 = groupsummary(Result_al, 'Date');
    datesToDelete = groupedTable2.Date(groupedTable2.GroupCount <= low_con);
    rowsToDelete = ismember(Result_al.Date, datesToDelete);
    Result_al(rowsToDelete, :) = [];

    elseif strcmpi(Third_filter, 'N')
    fprintf('Skip Filter 3 - Data Reliability \n');
end

save(fullfile(output, [output 'CS2_f3.mat']));

%% Making representative value by median

    %no_data = size(Result_al, 1);
    groupedTable3 = groupsummary(Result_al, 'Date', 'median', 'WL_CS2');
    %no_tracks = size(groupedTable3, 1);
    groupedTable4 = groupsummary(Result_al, 'Date', 'median', 'CI');
    
    %WL_Alt_al_med = groupedTable3.median_WL_CS2;
    %WL_Real_al_med = groupedTable4.median_gauged_WL;

%% Filter 4 : Applying Hydrological Filter

if strcmpi(Fourth_filter, 'Y')
    fprintf('Start Filter 4 - Hydrological Validity \n');

h = 1;
del_groupedTable3 = [];


while h < size(groupedTable3, 1) - 1
    
            fprintf('%d / %d completed', h, size(groupedTable3, 1) - 2)
    
    if groupedTable3.median_WL_CS2(h+1) > groupedTable3.median_WL_CS2(h)
        
        BefLv = groupedTable3.median_WL_CS2(h);
        PreLv = groupedTable3.median_WL_CS2(h+1);
        
        Shp_BefLv = shaperead(fullfile(shppath, sprintf(('%s_%s.shp'), DamName, num2str(round(BefLv)))));
        Shp_PreLv = shaperead(fullfile(shppath, sprintf(('%s_%s.shp'), DamName, num2str(round(PreLv)))));
        
        BefLv_X = Shp_BefLv.X;
        %BefLv_X = BefLv_X(~isnan(BefLv_X));
        BefLv_Y = Shp_BefLv.Y;
        %BefLv_Y = BefLv_Y(~isnan(BefLv_Y));
        [Bef_X_m, Bef_Y_m] = projfwd(utmCrs, BefLv_Y, BefLv_X);
        Bef_P = polyshape(Bef_X_m, Bef_Y_m, 'Simplify', false);
        Bef_Area_m2 = area(Bef_P);
        Bef_Area_km2 = Bef_Area_m2 / 1e6;
        
        PreLv_X = Shp_PreLv.X;
        %PreLv_X = PreLv_X(~isnan(PreLv_X));
        PreLv_Y = Shp_PreLv.Y;
        %PreLv_Y = PreLv_Y(~isnan(PreLv_Y));
        [Pre_X_m, Pre_Y_m] = projfwd(utmCrs, PreLv_Y, PreLv_X);
        Pre_P = polyshape(Pre_X_m, Pre_Y_m, 'Simplify', false);
        Pre_Area_m2 = area(Pre_P);
        Pre_Area_km2 = Pre_Area_m2 / 1e6;
        
        Delta_Sto = (Pre_Area_m2 + Bef_Area_m2)/2 * (PreLv - BefLv);
        
        %BefSto = dam_sto(BefLv, DamName);
        RF_stt = groupedTable3.Date(h) - days(3);
        RF_end = groupedTable3.Date(h+1) - days(1);
        idx = Gauged_daily.Gauged_time_daily >= RF_stt & Gauged_daily.Gauged_time_daily <= RF_end;
        RF = Gauged_daily(idx, :);
        Total_RF = sum(RF.Gauged_rainfall_daily);
        Inflow = Area * Total_RF * runoff_rate * 1e3;

            if Delta_Sto > Inflow
                Del_Row = groupedTable3(h+1, :);
                Del_Row.BefLv = BefLv;
                Del_Row.DeffLv = PreLv - BefLv;
                Del_Row.RF_stt = RF_stt;
                Del_Row.Total_RF = Total_RF;
                del_groupedTable3 = [del_groupedTable3; Del_Row];
                groupedTable3(h+1, :) = [];
            else
                h = h + 1;
            end
    else
        h = h + 1;
    end
    
    fprintf(repmat('\b', 1, numel(sprintf('%d / %d completed', h, size(groupedTable3, 1) - 1))));
    
end

groupedTable4 = groupedTable4(ismember(groupedTable4.Date, groupedTable3.Date), :);

    elseif strcmpi(Fourth_filter, 'N')
     fprintf('Skip Filter 4 - Hydrological Validity \n');

end

save(fullfile(output, [output 'CS2_f4.mat']));

%% The Result
    %fprintf('First_filter = %s, Second_filter = %s, Third_filter = %s, Fourth_filter = %s\n', First_filter, Second_filter, Third_filter, Fourth_filter);
    
    %Table_S_Dy1 = [];
    %no_cycle = 0;
    %no_data = 0;
            
    %no_data = sum(groupedTable3.GroupCount);
    %no_cycle = size(groupedTable3, 1);
    
    %idx = ismember(groupedTable4.Date, groupedTable3.Date);
    %groupedTable4 = groupedTable4(idx, :);
    
    %WL_Alt_al_med = groupedTable3.median_WL_CS2;
    %WL_Real_al_med = groupedTable4.median_gauged_WL;
            
    %mean_bias = mean(WL_Alt_al_med - WL_Real_al_med);
    %mean_bias = round(mean_bias, 4);
    %rmse = sqrt(mean((WL_Real_al_med - WL_Alt_al_med).^2));
    %rmse = round(rmse, 4);
    %ubrmse = sqrt(mean((WL_Real_al_med - WL_Alt_al_med).^2) - mean_bias.^2);
    %ubrmse = round(ubrmse, 4);
    %correlation = corr(WL_Real_al_med(~isnan(WL_Real_al_med) & ~isnan(WL_Alt_al_med)), WL_Alt_al_med(~isnan(WL_Real_al_med) & ~isnan(WL_Alt_al_med)));
    %correlation = round(correlation, 4);

    %mdl = fitlm(WL_Real_al_med, WL_Alt_al_med);
	%Ord_Rsq = mdl.Rsquared.Ordinary;
	%Ord_Rsq = round(Ord_Rsq, 4);
	%Adj_Rsq = mdl.Rsquared.Adjusted;
	%Adj_Rsq = round(Adj_Rsq, 4);

	%Stat_filename = ['!Stat_CS2' output '.csv'];
    %Dam = {DamName};
    %S_Dy1 = table(Dam, no_data, no_cycle, 'VariableNames', {'Dam', 'No_data', 'No_cycle'});
 	%S_Dy1 = table(CELL2TABLE(DamName), First_filter, Second_filter, Third_filter, no_data, no_cycle, ubrmse, Ord_Rsq, mean_bias, rmse, correlation, Adj_Rsq, 'VariableNames', {'DamName', 'First_filter', 'Second_filter', 'Third_filter', 'No_data', 'No_cycle', 'ubrmse', 'Ord_Rsq', 'mean_bias', 'rmse', 'correlation', 'Adj_Rsq'});
    %Table_S_Dy1 = [Table_S_Dy1; S_Dy1];
	%writetable(S_Dy1, Stat_filename);
    
    SatName = repmat({SatName}, height(groupedTable3), 1);
    Result_Graph = table(SatName, groupedTable3.Date, groupedTable3.GroupCount, groupedTable3.median_WL_CS2, groupedTable4.median_CI, ...
                     'VariableNames', {'SatName', 'Date', 'GroupCount', 'median_WL_CS2', 'median_CI'});
                 
    assignin('base', 'Result_CS2', Result_Graph);                 

    save(fullfile(output, [output 'CS2_ff.mat']));
    
%    clear Table_S_Dy1 del_groupedTable3 Result_al T2 DyWMs Sel_Set_FL Passfilename t_test* Sel_Set_* Close_Time Real_WL
    
fprintf('\nCryoSat-2 done\n\n')