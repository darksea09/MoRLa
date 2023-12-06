    fprintf('\n Copyright of MoRLa © UNSW and K-water 2023. All rights reserved. \n\n');
    
    % Input entering by user
    parameters.DamName = input('Please enter the dam name(e.g., CJ or DC): ', 's');
    parameters.Min_Lv = input('\nPlease enter the minimum operation level(m) of resorvoir and water mask(e.g., 111 or 61): ');
    parameters.Max_Lv = input('\nPlease enter the maximum operation level(m) of resorvoir and water mask(e.g., 144 or 79): ');
    parameters.path = input('\nPlease enter the directory path of satellite data: ', 's');
    parameters.shppath = input('\nPlease enter the directory path of water-mask data: ', 's');

    validInput = false;
    while ~validInput
        parameters.First_filter = input('\nWould you apply the first filter for Standard Operation Level?(Y/N): ', 's');
        if strcmpi(parameters.First_filter, 'Y') || strcmpi(parameters.First_filter, 'N')
            validInput = true;
        else
            warning('Invalid input for first filter. Please enter "Y" or "N" only.');
        end
    end

    validInput = false;
    while ~validInput
        parameters.Second_filter = input('\nWould you apply the second filter for Dynamic Water Mask? (Y/N): ', 's');
        if strcmpi(parameters.Second_filter, 'Y') || strcmpi(parameters.Second_filter, 'N')
            validInput = true;
        else
            warning('Invalid input for second filter. Please enter "Y" or "N" only.');
        end
    end
    
    validInput = false;
    while ~validInput
        parameters.Third_filter = input('\nWould you apply the third filter for Data Reliability? (Y/N): ', 's');
        if strcmpi(parameters.Third_filter, 'Y') || strcmpi(parameters.Third_filter, 'N')
            validInput = true;
        else
            warning('Invalid input for third filter. Please enter "Y" or "N" only.');
        end
    end
    
    if strcmpi(parameters.Third_filter, 'Y')
        validInput = false;
        while ~validInput
            parameters.low_con = input('Please enter the low-confidence value (1 or 2 or 3): ');
        if any(parameters.low_con == [1, 2, 3])
            validInput = true;
        else
            warning('Invalid input for low-confidence value. Please enter 1, 2, or 3 only.');
        end
        end
    else
        parameters.low_con = [];
    end
    
    validInput = false;
    while ~validInput
        parameters.Fourth_filter = input('\nWould you apply the fourth filter for Hydrological Validity? (Y/N): ', 's');
        if strcmpi(parameters.Fourth_filter, 'Y') || strcmpi(parameters.Fourth_filter, 'N')
            validInput = true;
        else
            warning('Invalid input for fourth filter. Please enter "Y" or "N" only.');
        end
    end
    
    if strcmpi(parameters.Fourth_filter, 'Y')
        validInput = false;
        while ~validInput
            parameters.Area = input('Please enter the basin area(km^2) (e.g., 6648, 3204): ');
            parameters.RGpath = input('Please enter the directory path of rainfall data: ', 's');
            parameters.runoff_rate = input('Please enter the runoff rate(0.2~0.9): ');
        if any(parameters.runoff_rate == [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
            validInput = true;
        else
            warning('Invalid input for runoff_rate. Please enter within the range 0.2~0.9 only.');
        end
        end
    else
        parameters.runoff_rate = [];
        parameters.Area = [];
    end
    
    parameters.numWorkers = input('\nPlease enter the number of cores you would like to use, considering your computer''s performance : ');
    parameters.output = input('\nPlease enter the directory path for the output: ', 's');
    
    % Fixed input
    parameters.Lv_range = parameters.Max_Lv:-1:parameters.Min_Lv;
    
% Call the converted functions with the defined parameters
S_AltiFunc_S3A(parameters.path, parameters.DamName, parameters.Lv_range, parameters.Area, ...
    parameters.First_filter, parameters.Second_filter, parameters.Third_filter, parameters.low_con, ...
    parameters.Fourth_filter, parameters.runoff_rate, parameters.output, parameters.RGpath, parameters.shppath, parameters.numWorkers);

S_AltiFunc_S3B(parameters.path, parameters.DamName, parameters.Lv_range, parameters.Area, ...
    parameters.First_filter, parameters.Second_filter, parameters.Third_filter, parameters.low_con, ...
    parameters.Fourth_filter, parameters.runoff_rate, parameters.output, parameters.RGpath, parameters.shppath, parameters.numWorkers);

S_AltiFunc_CS2(parameters.path, parameters.DamName, parameters.Lv_range, parameters.Area, ...
    parameters.First_filter, parameters.Second_filter, parameters.Third_filter, parameters.low_con, ...
    parameters.Fourth_filter, parameters.runoff_rate, parameters.output, parameters.RGpath, parameters.shppath, parameters.numWorkers);

S_AltiFunc_IS2(parameters.path, parameters.DamName, parameters.Lv_range, parameters.Area, ...
    parameters.First_filter, parameters.Second_filter, parameters.Third_filter, parameters.low_con, ...
    parameters.Fourth_filter, parameters.runoff_rate, parameters.output, parameters.RGpath, parameters.shppath, parameters.numWorkers);


Real_Gauged_daily = readtable(fullfile(parameters.RGpath, ['DailyRain_' parameters.DamName '.xlsx']));
Gauged_time_daily = datetime(Real_Gauged_daily.Var1,'TimeZone','UTC');
Gauged_rainfall_daily = Real_Gauged_daily.Var2;
Gauged_daily = table(Gauged_time_daily, Gauged_rainfall_daily);

% Loading each result from func.
Result_S3B.Properties.VariableNames = {'SatName', 'Date', 'GroupCount', 'WL_sat', 'CI'};
Result_S3A.Properties.VariableNames = {'SatName', 'Date', 'GroupCount', 'WL_sat', 'CI'};
Result_CS2.Properties.VariableNames = {'SatName', 'Date', 'GroupCount', 'WL_sat', 'CI'};
Result_IS2.Properties.VariableNames = {'SatName', 'Date', 'GroupCount', 'WL_sat', 'CI'};

All_Graph = table();
if istable(Result_S3B)
    All_Graph = [All_Graph; Result_S3B];
end
if istable(Result_S3A)
    All_Graph = [All_Graph; Result_S3A];
end
if istable(Result_CS2)
    All_Graph = [All_Graph; Result_CS2];
end
if istable(Result_IS2)
    All_Graph = [All_Graph; Result_IS2];
end

%All_Graph = [Result_S3B; Result_S3A; Result_CS2; Result_IS2];
All_Graph.Date = datetime(All_Graph.Date, 'InputFormat', 'dd-MMM-yyyy');
All_Graph = sortrows(All_Graph, 'Date');

years = year(All_Graph.Date);
months = month(All_Graph.Date);

month_strs = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

old_year = min(years);
old_month = min(months(years == old_year));
old_month = month_strs{old_month};

recent_year = max(years);
recent_month = max(months(years == recent_year));
end_day = eomday(recent_year, recent_month);
recent_month = month_strs{recent_month};


% 월별 강우량 계산
Gauged_daily.Month = month(Gauged_daily.Gauged_time_daily);
Gauged_daily.Year = year(Gauged_daily.Gauged_time_daily);
monthlyRainAll = groupsummary(Gauged_daily, {'Year', 'Month'}, 'sum', 'Gauged_rainfall_daily');

% 그래프 생성
figure('Position', [100, 100, 1260, 600]);

yyaxis left;

plot(datenum(All_Graph.Date), All_Graph.WL_sat, '-b', 'LineWidth', 1.0, 'HandleVisibility', 'off');
hold on;
S3B_Data = All_Graph(strcmp(All_Graph.SatName, 'S3B'), :);
S3A_Data = All_Graph(strcmp(All_Graph.SatName, 'S3A'), :);
CS2_Data = All_Graph(strcmp(All_Graph.SatName, 'CS2'), :);
IS2_Data = All_Graph(strcmp(All_Graph.SatName, 'IS2'), :);

if ~isempty(S3B_Data)
    errorbar(datenum(S3B_Data.Date), S3B_Data.WL_sat, S3B_Data.CI, 'o', 'Color', 'black', 'MarkerFaceColor', [0.8, 0.6, 0.7], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'S3B');
end

if ~isempty(S3A_Data)
    errorbar(datenum(S3A_Data.Date), S3A_Data.WL_sat, S3A_Data.CI, 'o', 'Color', 'black', 'MarkerFaceColor', [0.6, 0.8, 0.7], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'S3A');
end

if ~isempty(CS2_Data)
    errorbar(datenum(CS2_Data.Date), CS2_Data.WL_sat, CS2_Data.CI, 'o', 'Color', 'black', 'MarkerFaceColor', [0.7, 0.8, 0.9], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'CS2');
end

if ~isempty(IS2_Data)
    errorbar(datenum(IS2_Data.Date), IS2_Data.WL_sat, IS2_Data.CI, 'o', 'Color', 'black', 'MarkerFaceColor', [0.9, 0.7, 0.8], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'IS2');
end

% x-axis
datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
ylabel('Water Level (m)', 'FontWeight', 'bold', 'FontSize', 14);
set(gca, 'YLim', [parameters.Min_Lv, parameters.Max_Lv+10], 'YTick', parameters.Min_Lv:5:parameters.Max_Lv, 'FontWeight', 'bold', 'FontSize', 14);

% y-axis for right
yyaxis right;
monthlyRain = monthlyRainAll.sum_Gauged_rainfall_daily;
monthlyDates = datetime(monthlyRainAll.Year, monthlyRainAll.Month, 15); % 월 중간을 대표 날짜로 사용
maxRainValue = max(monthlyRain);

% rainfall stick
bar(datenum(monthlyDates), monthlyRain, 'FaceColor', [0.8, 0.8, 1.0], 'DisplayName', 'Monthly Rainfall');

% labelling each month
for i = 1:length(monthlyDates)
    if year(monthlyDates(i)) == 2022
        text(datenum(monthlyDates(i)), monthlyRain(i)+10, datestr(monthlyDates(i), 'mmm'), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'Rotation', 20, 'FontWeight', 'bold', 'FontSize', 14);
    end
end

% Attribution of graph
set(gca, 'YDir', 'reverse', 'YLim', [0, 2*maxRainValue+300], 'YTick', 0:150:maxRainValue+100);
ylabel('Rainfall (mm)', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'southeast', 'FontWeight', 'bold', 'FontSize', 14);
grid off;
xlim([datenum(sprintf('01-%s-%d', old_month, old_year)), datenum(sprintf('%02d-%s-%d', end_day, recent_month, recent_year))]);

saveas(gcf, fullfile(parameters.output, sprintf('%s_WLv_Time_Series.jpg', parameters.output)));

All_Graph.Monthly_Rainfall = NaN(size(All_Graph, 1), 1);
for i = 1:height(All_Graph)
    currentDate = datetime(All_Graph.Date(i), 'InputFormat', 'dd-MMM-yyyy');
    current_year = year(currentDate);
    current_month = month(currentDate);
    
    matched_row = (monthlyRainAll.Year == current_year) & (monthlyRainAll.Month == current_month);
    
    if any(matched_row)
        All_Graph.Monthly_Rainfall(i) = monthlyRainAll{matched_row, 4};
    end
end

writetable(All_Graph, fullfile(parameters.output, ['Final_' parameters.output '.csv']));

save(fullfile(parameters.output, [parameters.output 'Final.mat']));