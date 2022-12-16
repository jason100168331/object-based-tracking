
clc
clear all
for yr = 51:70
    year = int2str(1949+yr);
    for mon=1:12
        Inpath = strcat('/Volumes/IMERG/ERE/ERA5/',year,'/',year,'_',int2str(mon),'.mat');
        Outpath = strcat('/Volumes/IMERG/ERE/ERA5/',year,'/',year,'_',int2str(mon),'_merged_.mat');
        load(Inpath)
        
        extreme_event_west_boundary = [];
        extreme_event_east_boundary = [];
        boundary_id = [];
        min_lon = [];
        max_lon = [];
        
        east_bound = 178.75;
        west_bound = -178.75;

        for i=1:length(extreme_event_all)
            pos=i;
            lon = [extreme_event_all(pos).Lon];
            lat = [extreme_event_all(pos).Lat];
            precip = [extreme_event_all(pos).Precip];
            lon = double(lon);
            lat = double(lat);
            lon = lon/100;
            lat = lat/100;
            

            if max(lon) >= east_bound 
                extreme_event_east_boundary=[extreme_event_east_boundary,extreme_event_all(pos)];
                boundary_id = [boundary_id,pos];
            elseif  min(lon)<= west_bound
                extreme_event_west_boundary=[extreme_event_west_boundary,extreme_event_all(pos)];
                boundary_id = [boundary_id,pos];
            else
                %boundary_id = [boundary_id,0];
            end

            min_lon = [min_lon,min(lon)];
            max_lon = [max_lon,max(lon)];
            precip = precip./100;
            %lon(precip<3)=[];
            %lat(precip<3)=[];
            time = [extreme_event_all(pos).startSecond]';
            startTime = datestr(datetime(min(time), 'ConvertFrom', 'posixtime'));
            endTime = datestr(datetime(max(time), 'ConvertFrom', 'posixtime'));                        
        end

        extreme_event_all(boundary_id)=[];
       
        poss = 0;
        extreme_event_boundary = struct('id',[],'Lat',[],'Lon',[],'Precip',[],'startSecond',[]);
        merged_id = [];
        east_merged_id = [];
        west_merged_id = [];

        for i=1:length(extreme_event_west_boundary)
            lon_west = [extreme_event_west_boundary(i).Lon];
            lat_west = [extreme_event_west_boundary(i).Lat];
            precip_west = [extreme_event_west_boundary(i).Precip];
            startSecond_west = [extreme_event_west_boundary(i).startSecond];

            lon_west = double(lon_west);
            lat_west = double(lat_west);
            lon_west = lon_west./100;
            lat_west = lat_west./100;

            lat_west(lon_west>west_bound)=[];
            startSecond_west(lon_west>west_bound)=[];
            lon_west(lon_west>west_bound)=[];

            for j=1:length(extreme_event_east_boundary)
                lon_east = [extreme_event_east_boundary(j).Lon];
                lat_east = [extreme_event_east_boundary(j).Lat];
                precip_east = [extreme_event_east_boundary(j).Precip];
                startSecond_east = [extreme_event_east_boundary(j).startSecond];

                lon_east = double(lon_east);
                lat_east = double(lat_east);
                lon_east = lon_east./100;
                lat_east = lat_east./100;

                lat_east(lon_east<east_bound)=[];
                startSecond_east(lon_east<east_bound)=[];
                lon_east(lon_east<east_bound)=[];

                if sum(ismember(lat_east,lat_west))>0 & sum(ismember(startSecond_east,startSecond_west))>0
                    poss = poss+1;
                    id_temp = [i,j];
                    merged_id = [merged_id;id_temp];
                    west_merged_id = [west_merged_id,i];
                    east_merged_id = [east_merged_id,j];                                
                end
            end
        end


        east_merged_id_uni = unique(east_merged_id);
        west_merged_id_uni = unique(west_merged_id);
        west_set = cell(length(west_merged_id_uni),2);

        for i=1:length(west_merged_id_uni)
            temp = find(west_merged_id==west_merged_id_uni(i));
            west_set(i,1) = {west_merged_id_uni(i)};
            west_set(i,2) = {east_merged_id(temp)};
        end

        test=cell(1,1);
        posss=0;
        for i=1:size(merged_id,1)
            temp = find(east_merged_id==east_merged_id(i));
            if length(temp) > 1
                posss=posss+1;
                test{posss,1} = west_merged_id(temp);
            end
        end

        east_id_have_merged=[];

        for i=1:size(west_set,1)
            lon_west = [extreme_event_west_boundary(west_set{i,1}).Lon];
            lat_west = [extreme_event_west_boundary(west_set{i,1}).Lat];
            precip_west = [extreme_event_west_boundary(west_set{i,1}).Precip];
            startSecond_west = [extreme_event_west_boundary(west_set{i,1}).startSecond];

            lon_merged = lon_west;
            lat_merged = lat_west;
            precip_merged = precip_west;
            second_merged = startSecond_west;

            target_id = west_set{i,2};

            for j=1:length(target_id)
                if sum(ismember(target_id(j),east_id_have_merged))==0
                    lon_east = [extreme_event_east_boundary(target_id(j)).Lon];
                    lat_east = [extreme_event_east_boundary(target_id(j)).Lat];
                    precip_east = [extreme_event_east_boundary(target_id(j)).Precip];
                    startSecond_east = [extreme_event_east_boundary(target_id(j)).startSecond];

                    lon_merged = [lon_merged;lon_east];
                    lat_merged = [lat_merged;lat_east];
                    precip_merged = [precip_merged;precip_east];
                    second_merged = [second_merged;startSecond_east];

                    east_id_have_merged = [east_id_have_merged,target_id(j)];
                end
            end

            extreme_event_west_boundary(west_set{i,1}).Lon = lon_merged;
            extreme_event_west_boundary(west_set{i,1}).Lat = lat_merged;
            extreme_event_west_boundary(west_set{i,1}).Precip = precip_merged;
            extreme_event_west_boundary(west_set{i,1}).startSecond = second_merged;

        end

        out = unique(cellfun(@num2str,test,'uni',0));
        west_delete_id=[];
        if length(test{1,1})>0
            for i=1:length(out)
                west_target_id = out{i};
                west_target_id = strsplit(west_target_id,'  ');

                lon_merged = [];
                lat_merged = [];
                precip_merged = [];
                second_merged = [];

                for j=1:length(west_target_id)
                    target_id = str2double(west_target_id{j});

                    lon_west = [extreme_event_west_boundary(target_id).Lon];
                    lat_west = [extreme_event_west_boundary(target_id).Lat];
                    precip_west = [extreme_event_west_boundary(target_id).Precip];
                    startSecond_west = [extreme_event_west_boundary(target_id).startSecond];

                    if j>1
                        west_delete_id = [west_delete_id,target_id];
                    end
                end
            end
        end

        extreme_event_west_boundary(west_delete_id)=[];
        extreme_event_east_boundary(east_id_have_merged)=[];

        for i=1:length(extreme_event_west_boundary)
            fprintf('i= %d\n',i);
            lon = [extreme_event_west_boundary(i).Lon];
            lat = [extreme_event_west_boundary(i).Lat];
            precip = [extreme_event_west_boundary(i).Precip];
            startTime_second_event = [extreme_event_west_boundary(i).startSecond]';

            lon_event = double(lon);
            lat_event = double(lat);
            precip = double(precip);
            lon_event = lon_event./100;
            lat_event = lat_event./100;
            precip_event = precip./100;
            
            lon_360 = lon_event;
            
            for j=1:length(lon_360)
                if lon_360(j)<0
                    lon_360(j) = lon_360(j)+360;
                end
            end

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

            for j=1:length(unique_startTime)
                num = find(startTime_second_event == unique_startTime(j));
                lat_temp = mean(lat_event(num));
                lon_temp = mean(lon_360(num));
                centroid_lat_list(j) = lat_temp/100;
                centroid_lon_list(j) = lon_temp/100;
                instant_area(j) = length(num);
            end

            storm_track = [centroid_lat_list,centroid_lon_list];

            distance_all = 0;
            for ii=2:length(centroid_lat_list)
                distance = distt(centroid_lat_list(ii-1),centroid_lon_list(ii-1),centroid_lat_list(ii),centroid_lon_list(ii));
                distance_all = distance_all + distance;
            end

            speed = distance_all/duration;
        
            volume = length(lat_event);
            area = length(unique_lat);  

            extreme_event_west_boundary(i).id = i;
            extreme_event_west_boundary(i).Lat = int16(lat_event.*100);
            extreme_event_west_boundary(i).Lon = int16(lon_event.*100);
            extreme_event_west_boundary(i).Precip = precip_event_100;
            extreme_event_west_boundary(i).startSecond = startTime_second_event; 
            extreme_event_west_boundary(i).Area = area;
            extreme_event_west_boundary(i).Volume = volume;
            extreme_event_west_boundary(i).instantMaxArea = max(instant_area);
            extreme_event_west_boundary(i).Duration = duration;
            extreme_event_west_boundary(i).Speed = speed;            
            extreme_event_west_boundary(i).maxPrecip = maxPrecip;           
            extreme_event_west_boundary(i).totalPrecip = totalPrecip;
            extreme_event_west_boundary(i).meanPrecip = meanPrecip;   
        end

        extreme_event_all = [extreme_event_all,extreme_event_west_boundary,extreme_event_east_boundary];
        save(Outpath,'extreme_event_all');
        
    end
end