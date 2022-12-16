clc
clear all

for yr = 2:70
    year = 1978+yr;
    fprintf('year= %d\n',year);
        
    if mod(year,4) ~= 0
        range = cell(12,1);
        range{1} = [1:744];range{2} = [745:1416];
        range{3} = [1417:2160];range{4} = [2161:2880];
        range{5} = [2881:3624];range{6} = [3625:4344];
        range{7} = [4345:5088];range{8} = [5089:5832];
        range{9} = [5833:6552];range{10} = [6553:7296];
        range{11} = [7297:8016];range{12} = [8017:8760];
        
    elseif mod(year,4) == 0
        range = cell(12,1);
        range{1} = [1:744];range{2} = [745:1440];
        range{3} = [1441:2184];range{4} = [2185:2904];
        range{5} = [2905:3648];range{6} = [3649:4368];
        range{7} = [4369:5112];range{8} = [5113:5846];
        range{9} = [5847:6762];range{10} = [6763:7320];
        range{11} = [7321:8040];range{12} = [8041:8784];
    end

    Inpath=strcat('/Volumes/IMERG/ERA5/Total Precipitation/',int2str(year),'/download.nc');
    %Indir=dir(fullfile(Inpath,'*.HDF'));
    %Inlength=length(Indir);

    threshold = 1;
    event_size_threshold = 50;
    radii = 0.75;
    filter = ones(int32(radii/0.25),int32(radii/0.25));

    lat_all = ncread(Inpath,'latitude');
    lon_all = ncread(Inpath,'longitude');
    lon_all = lon_all-180;

    lat_all_left = linspace(-90-(size(filter,1)-1)*0.25,-90,size(filter,1));
    lat_all_left(1)=[];
    lat_all_right = linspace(90,90+(size(filter,1)-1)*0.25,size(filter,1));
    lat_all_right(end)=[];
    lat_all = [lat_all_left';lat_all;lat_all_right'];


    lon_all_left = linspace(-180-(size(filter,1)-1)*0.25,-180,size(filter,1));
    lon_all_left(1)=[];
    lon_all_right = linspace(180,180+(size(filter,1)-1)*0.25,size(filter,1));
    lon_all_right(end)=[];
    lon_all = [lon_all_left';lon_all;lon_all_right'];

    lat_convolved = [];
    for i=1:(length(lat_all)-(size(filter,1)-1))
        lat_sum = sum(lat_all(i:i+(size(filter,1)-1)));
        lat_mean = lat_sum/size(filter,1);
        lat_convolved = [lat_convolved;lat_mean];
    end

    lon_convolved = [];
    for i=1:(length(lon_all)-(size(filter,1)-1))
        lon_sum = sum(lon_all(i:i+(size(filter,1)-1)));
        lon_mean = lon_sum/size(filter,1);
        lon_convolved = [lon_convolved;lon_mean];
    end

    for mon = 1:12
        fprintf('mon= %d\n',mon);
        ind = 0;
        
        time = ncread(Inpath,'time',range{mon}(1),length(range{mon}));
        
        startTime_second_all = [];
        
        for i=1:length(time)

            fprintf('i= %d\n',i);
            format_out = 'yyyy-mm-dd HH:MM:SS';
            startTime_second = posixtime(datetime(datestr(datenum(double(time(i))/24)+datenum(1900,01,01),format_out)));
            %startTime_second = datetime(datestr(datenum(double(time(i))/24)+datenum(1900,01,01)));

            startTime_second_all = [startTime_second_all,startTime_second];
        end

        extreme_inventory_uncal = zeros(1442,723,length(time));
        precip_all = ncread(Inpath,'tp',[1,1,range{mon}(1)],[1440,721,length(time)]);

        for p=1:length(time)
            ind = ind + 1;
            fprintf('ind= %d\n',ind);
            precip_id  = range{mon}(p);
            fprintf('precip_id= %d\n',precip_id);
            precip = precip_all(:,:,p);
            precip = [precip(721:1440,:);precip(1:720,:)];
            
            precip = precip.*1000;
            precip(precip<0)=0;

            precip_convolved = conv2(precip,filter,'full');
            precip_convolved = precip_convolved./length(filter(:));

            precip_convolved(precip_convolved<threshold)=0;                                               

            extreme_inventory_uncal(:,:,ind) = precip_convolved;
        end

        connectiveEvent = bwlabeln(extreme_inventory_uncal,6); 

        %clearvars extreme_inventory_uncal;

        pixelList = regionprops(connectiveEvent,'PixelList');
        
        f=0;

        extreme_event_all = struct('id',[],'Lat',[],'Lon',[],'Precip',[],'startSecond',[],...
            'Volume',[],'Area',[],'instantMaxArea',[],'Duration',[],'Speed',[],...
            'maxPrecip',[],'totalPrecip',[],'meanPrecip',[]);

        for i1=1:length(pixelList)        

            fprintf('i1 = %d\n',i1);
            pixels = pixelList(i1).PixelList;
            if length(pixels)>event_size_threshold
                f=f+1;

                x = pixels(:,1);
                y = pixels(:,2);
                z = pixels(:,3);

                lat_event = zeros(size(pixels,1),1);
                lon_event = zeros(size(pixels,1),1);
                precip_event = zeros(size(pixels,1),1);

                startTime_second_event = zeros(size(pixels,1),1);

                for i2=1:size(pixels,1)
                    lat_event(i2) = lat_convolved(x(i2));
                    lon_event(i2) = lon_convolved(y(i2));
                    precip_event(i2) = extreme_inventory_uncal(y(i2),x(i2),z(i2));
                    startTime_second_event(i2) = startTime_second_all(z(i2));
                end

                lat_event = lat_event.*100;
                lon_event = lon_event.*100;
                unique_cor = unique([double(lat_event');double(lon_event')]','rows');
                unique_lat = unique_cor(:,1);
                unique_lon = unique_cor(:,2);
                precip_event = round(precip_event,2);
                precip_event_100 = precip_event.*100;
                precip_event_100 = round(precip_event_100);
                maxPrecip = max(precip_event_100);
                totalPrecip = sum(precip_event_100)/2;
                meanPrecip = mean(precip_event_100);

                lat_event = int16(lat_event);
                lon_event = int16(lon_event);
                precip_event_100 = int16(precip_event_100);

                duration = max(startTime_second_event) - min(startTime_second_event);

                if duration == 0
                    duration = 1;
                else
                    duration = (duration/3600) + 1;
                end

                unique_startTime = unique(startTime_second_event);

                centroid_lat_list = zeros(length(unique_startTime),1);
                centroid_lon_list = zeros(length(unique_startTime),1);
                instant_area = zeros(length(unique_startTime),1);

                for i=1:length(unique_startTime)
                    num = find(startTime_second_event == unique_startTime(i));
                    lat_temp = mean(lat_event(num));
                    lon_temp = mean(lon_event(num));
                    centroid_lat_list(i) = lat_temp/100;
                    centroid_lon_list(i) = lon_temp/100;
                    instant_area(i) = length(num);
                end

                storm_track = [centroid_lat_list,centroid_lon_list];
                
                % speed
                distance_all = 0;
                for i=2:length(centroid_lat_list)
                    distance = distt(centroid_lat_list(i-1),centroid_lon_list(i-1),centroid_lat_list(i),centroid_lon_list(i));
                    distance_all = distance_all + distance;
                end
                
                speed = distance_all/duration;

                volume = length(lat_event);
                area = length(unique_lat);  

                extreme_event_all(f).id = f;
                extreme_event_all(f).Lat = lat_event;
                extreme_event_all(f).Lon = lon_event;
                extreme_event_all(f).Precip = precip_event_100;
                extreme_event_all(f).startSecond = startTime_second_event; 
                extreme_event_all(f).Area = area;
                extreme_event_all(f).Volume = volume;
                extreme_event_all(f).instantMaxArea = max(instant_area);
                extreme_event_all(f).Duration = duration;
                extreme_event_all(f).Speed = speed;            
                extreme_event_all(f).maxPrecip = maxPrecip;           
                extreme_event_all(f).totalPrecip = totalPrecip;
                extreme_event_all(f).meanPrecip = meanPrecip;        
            end
        end 

        volume_all = [extreme_event_all.Volume];
        [volume_all,sortIdx] = sort(volume_all,'descend');
        extreme_event_all = extreme_event_all(sortIdx);

        Outpath = strcat('/Volumes/IMERG/ERE/ERA5/',int2str(year),'/',...
            int2str(year),'_',int2str(mon),'.mat');
        save(Outpath,'extreme_event_all');        
    end
end