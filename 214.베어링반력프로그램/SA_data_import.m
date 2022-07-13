clear all; close all; clc;

%% SA data import 
[file, path] = uigetfile({'*.xlsx'});
file_ext = fullfile(path,file);

[~, ~, SA_data] = xlsread(file_ext,'Sheet1');
SA_data(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),SA_data)) = {''};

%% Shaft data extraction

for i=1:length(SA_data)
    Shaft_data_st(i)  = regexp(SA_data(i),'Shaft elements');
    Shaft_data_end(i) = regexp(SA_data(i),'Total');
end

Shaft_st_no  = find(cellfun('isempty', Shaft_data_st ) == 0 );
Shaft_end_no = find(cellfun('isempty', Shaft_data_end) == 0 );

Shaft_data=SA_data(Shaft_st_no+3:Shaft_end_no,1);
Shaft_data_no=regexp(Shaft_data,['\d+\.?\d*'],'match');

% ------ find max element size due to no. in the shaft element name 
for i=1:length(Shaft_data_no)
    temp_no(i) = length ( str2double(Shaft_data_no{i}) );
end
temp_no_max = max(temp_no);

% ------ reorder the shaft element starting from the back(max element no)
for i=1:length(Shaft_data_no)
    Shaft_element_size = length ( str2double(Shaft_data_no{i}) );
    if max(temp_no_max) > Shaft_element_size
        st_no = temp_no_max - Shaft_element_size + 1;
        Shaft_element(i, st_no : temp_no_max) = str2double(Shaft_data_no{i});
    else
        Shaft_element(i,:) = str2double(Shaft_data_no{i});
    end   
end

%% find R,BNNT,Prop,N,AFT

L_Aft_raw=strfind(Shaft_data,'Aft seal');
L_Aft_raw_log=~cellfun('isempty',L_Aft_raw);
L_Aft_idx=find(L_Aft_raw_log);
L_Aft=Shaft_element(L_Aft_idx(:,1,1),2);
L_Aft=mean(Shaft_element(L_Aft_idx(:,1,1),2));
%%Shaft_data_no=regexp(Shaft_data,['\w+\.?\d*'],'match');


%% Prop extraction

for i=1:length(SA_data)
    prop_st(i)  = regexp(SA_data(i),'PROPELLER - Non condition dependent');
    prop_end(i) = regexp(SA_data(i),'PROPELLER - Condition dependent');
end

prop_st_no  = find(cellfun('isempty', prop_st ) == 0 );
prop_end_no = find(cellfun('isempty', prop_end) == 0 );

prop_data=SA_data(prop_st_no+3:prop_end_no,1);
prop_data_no=regexp(prop_data,['\d+\.?\d*'],'match');

for i=1:length(prop_data_no)
    temp_no(i) = length ( str2double(prop_data_no{i}) );
end
temp_no_max = max(temp_no);

for i=1:length(prop_data_no)
    prop_size = length ( str2double(prop_data_no{i}) );
    if max(temp_no_max) > prop_size
        st_no = temp_no_max - prop_size + 1;
        prop(i, st_no : temp_no_max) = str2double(prop_data_no{i});
    else
        prop(i,:) = str2double(prop_data_no{i});
    end   
end

prop([find(prop==0)])=[];
L_prop=prop(1,1);
W_prop=prop(1,2);


%% Other case
% teststr = 'Test setup: MaxDistance = 60 m, Rate = 1.000, Permitted Error = 50 Operator Note:  Air Temperature=20 C, Wind Speed 16.375m/s, Altitude 5km (Cloudy)';
% regexp(teststr,[\d])
% regexp(teststr,['\d'])
% regexp(teststr,['\d'],'match')
% regexp(teststr,['\d+'],'match')
% regexp(teststr,['\d+.?'],'match')
% regexp(teststr,['\d+\.?'],'match')
% regexp(teststr,['\d+\.?\d?'],'match')
% regexp(teststr,['\d+\.?\d+?'],'match')
% regexp(teststr,['\d+\.?\d*?'],'match')
% regexp(teststr,['\d+\.?\d?'],'match')
% regexp(teststr,['\d+\.?\d*'],'match')
% 
% t2 = cell2mat(cellfun(@str2num,regexp(T, '(\d+) (\d+) (\d+) (\d+) (\d+) ', 'tokens','once'),'uniform',0))
% 
% 