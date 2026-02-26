% Julia Araujo - AUG2023
% Regridding glider data
% Addapted from script provided by Amandine Schaeffer and Saranya Kumar,
% with the help of Daneeja Mawren

% addpath(genpath('/home/jdearaujo/matlab_runs/my_packages'))

clear
clc

%% Parameters
% Path of the .csv file with the list of gliders that Amandine made
path_csv='C:\Users\Julia\Documents\projects\mhw_australia\data\gliders_database_final.csv';
path_nc= 'C:\Users\Julia\Documents\projects\mhw_australia\data\gliders_qc-regrid_vf\';

% Variables name you want to extract
dims_names={'LONGITUDE','LATITUDE','TIME','DEPTH'};
vars_names={'TEMP','PSAL','PRES','DOX2','CPHL'};

% Variables attributes
vars_attrb={'Temperature','Salinity','Pressure','Dissolved oxygen','Chlorophyll'};
vars_units={'Celsius','none','dbar','umol kg-1','mg m-3'};

% OPeNDAR urls (don't need to change this part)
url_catalog='https://thredds.aodn.org.au/thredds/catalog/IMOS/ANFOG/';
url_data='https://thredds.aodn.org.au/thredds/dodsC/IMOS/ANFOG/';

% Set depth resolution (every 1m depth)
depth_res=1;
max_depth=150; % Maximum depth of the dataset

% Set distance resolution (every 1km)
dist_res=1;

% Set bathymetry range for QC2
bathy_tasnsw=[20 120];
bathy_qldwas=[20 80];

% Set bottom boundary layer for QC3
bbl_margin=20;
shallow_glider_tasnsw=40;
shallow_glider_qldwas=30;

% Set window size for QC4 in chlorophyll
window=1000;

%% Load bathymetry 
bathy_lon=double(ncread('data/bathy_dbdb2_v30.nc','lon'));
bathy_lat=double(ncread('data/bathy_dbdb2_v30.nc','lat'));
bathy_height=double(ncread('data/bathy_dbdb2_v30.nc','height'));

%% Get the glider names (folder names of OPeNDAP link)
list_csv=readtable(path_csv);
names=list_csv.Name;
regions=list_csv.Region;
types=list_csv.ID;
types=cellfun(@(v)v(1:2),types,'UniformOutput',false);

%% Loop for every glider name (folder)
pb=waitbar(0,'Starting');

for i=117:length(names)
    disp(' ')
    disp(['File ' names{i}])
    disp(' ')
    clear DATA

    if ~exist([path_nc names{i} '_QCs.nc'],'file') % check if file already exists

        % Define thresholds for bathymetry range (QC2), BBL (QC3) and shallow gliders
        % for NSW/TAS and QLD/WAS
        if ismember(regions{i},{'TAS','NSW'})
            bathy_range=bathy_tasnsw;
            shallow_glider=shallow_glider_tasnsw;
        elseif ismember(regions{i},{'QLD','SW WA'})
            bathy_range=bathy_qldwas;
            shallow_glider=shallow_glider_qldwas;
        end
    
        % Select the glider type to follow the right OPeNDAP folder
        switch types{i}
            case 'SL'
                glider_type='slocum_glider/';
            case 'SG'
                glider_type='seaglider/';
                disp(['Seaglider: ' names{i}])
        end
    
        % List the files inside the OPeNDAP folder
        nc_filename=ls_opendap([url_catalog glider_type names{i} '/catalog.html']);
    
        %% Load dimensions (longitude, latitude, depth and time)
        disp('Loading dimensions...')
        for ii=1:length(dims_names)
            dims.(dims_names{ii})=ncread([url_data glider_type names{i} '/' nc_filename{1}],dims_names{ii});
        end
    
        %% Find peaks of depths (down-up casts)
        disp('Finding peaks...')
        [~,ind_pd]=findpeaks(dims.DEPTH,'MinPeakProminence',5);
        [~,ind_pu]=findpeaks(-dims.DEPTH,'MinPeakProminence',5);
    
        % Get first and last 'up peak'
        dummy_i=find(~isnan(dims.DEPTH),1,'first');
        dummy_f=find(~isnan(dims.DEPTH),1,'last');
        ind_pu=cat(1,dummy_i,ind_pu,dummy_f);
    
        % Sometimes, the function 'findpeaks' doesn't account for the right
        % 5-prominence value. Below, I correct this case (two cases found so far)
    
        % Get the wrong index
        if length(ind_pu)==length(ind_pd)
            ind_diff=ind_pd-ind_pu;
            wrong_ind=find(ind_diff<0,1,'first');
            ind_pd(wrong_ind)=[];
        end
    
        % Sometimes the function 'findpeaks' doesn't get the minimum peaks. 
        % Bellow, I correct this case (one case found so far)
    
        % Use the 'islocalmin' and 'islocalmax' with a lower prominence value
        if length(ind_pu)~=length(ind_pd)+1
            ind_pd=islocalmax(dims.DEPTH,'MinProminence',1);
            ind_pu=islocalmin(dims.DEPTH,'MinProminence',1);
            ind_pd=find(ind_pd);
            ind_pu=find(ind_pu);
    
            % Get first and last 'up peak'
            dummy_i=find(~isnan(dims.DEPTH),1,'first');
            dummy_f=find(~isnan(dims.DEPTH),1,'last');
            ind_pu=cat(1,dummy_i,ind_pu,dummy_f);
        end
    
        clear dummy_i dummy_f
    
        %% Binning dimensions to fixed depths
        disp('Binning dimensions to fixed depths...')
        cast=[];            
    
        for ii=1:length(dims_names)
            cast=casts(cast,dims.(dims_names{ii}),dims_names{ii},dims.DEPTH, ...
                       'down',ind_pu,ind_pd,max_depth,depth_res);
            cast=casts(cast,dims.(dims_names{ii}),dims_names{ii},dims.DEPTH, ...
                       'up',ind_pu,ind_pd,max_depth,depth_res);
        end
    
        %% Average down-up casts for dimensions
        disp('Averaging down-up casts for dimensions...')
        for ii=1:length(dims_names)
            var_dummy=cat(3,cast.down.(dims_names{ii}),cast.up.(dims_names{ii}));
            cast.new.(dims_names{ii})=mean(var_dummy,[1 3],'omitnan');        
        end
    
        cast.new=rmfield(cast.new,'DEPTH');
    
        var_dummy=cat(3,cast.down.DEPTH_ref,cast.up.DEPTH_ref);
        cast.new.DEPTH=mean(var_dummy,[2 3],'omitnan');
    
        disp(size(var_dummy));
    
        clear var_dummy

        %% Fill gaps in coordinates
        cast.new.LONGITUDE=fillgaps(cast.new.LONGITUDE);
        cast.new.LATITUDE=fillgaps(cast.new.LATITUDE);

        %% Clear empty coordinates
        clean_coords=mean([cast.new.LONGITUDE; cast.new.LATITUDE],1);
        is_nan=isnan(clean_coords);
        fields=fieldnames(cast.new);
        for ii=1:length(fields)-1
            f=fields{ii};
            cast.new.(f)(is_nan)=[];
        end
    
        %% Binning dimensions to fixed distances
        % Calculate distance in km
        dist=sw_dist(cast.new.LATITUDE,cast.new.LONGITUDE,'km');
        cast.new.DISTANCE=[0 cumsum(dist)];
        max_dist=ceil(max(cast.new.DISTANCE));
    
        disp('Binning dimensions to fixed distances...')
        for ii=1:length(dims_names)
            if strcmp(dims_names{ii},'DEPTH')
                DATA.(dims_names{ii})=cast.new.(dims_names{ii});
            else
                [new_distance,DATA.(dims_names{ii})]=binning(cast.new.DISTANCE, ...
                                                             cast.new.(dims_names{ii}), ...
                                                             max_dist,dist_res);
                DATA.DISTANCE=new_distance';
                DATA.(dims_names{ii})=DATA.(dims_names{ii})';
            end
        end
        
        %% Bathymetry flags (QC2)
        % Getting only depths from 20 to 80 or 120 meters
        bathy_data=[];
        
        % Find bathymetry for each point
        for ip=1:length(dims.LATITUDE)
            ind_lat=find(bathy_lat>dims.LATITUDE(ip),1,'first');
            ind_lon=find(bathy_lon>dims.LONGITUDE(ip),1,'first');
    
            if ~isempty(ind_lon)
                bathy_data(ip)=bathy_height(ind_lon,ind_lat);
            else
                bathy_data(ip)=NaN;
            end
        end
    
        bathy_data=abs(bathy_data);
        flags_bathy=(bathy_data>=bathy_range(1) & bathy_data<=bathy_range(2));
        
        %% BBL flags (bottom boundary layer) (QC3)
        glider_depth=dims.DEPTH;
        flags_bbl=glider_depth>(bathy_data-bbl_margin)';
    
        %% Quenching flags (QC5)
        % Getting only the daytime data
        time_vec=double(dims.TIME)+datetime(1950,1,1); % Add epoch offset
        time_local=time_vec+hours(10);                 % This is for reference
        
        % Vectorized solar position calculation
        [solar_az,solar_el]=SolarAzEl(time_vec,dims.LATITUDE,dims.LONGITUDE,0);
        
        % Logical indexing for daytime
        flags_daytime=solar_el>0;

        %% Loop for every variable
        for ii=1:length(vars_names)
            disp(' ')
            clear data
    
            % Sometimes, the variable doesn't exist in the file
            % In this case, the code will skip the variable
    
            try
                %% Load data
                disp(['Loading variable ' vars_names{ii} '...'])
                data.(vars_names{ii})=ncread([url_data glider_type names{i} '/' nc_filename{1}],vars_names{ii});
                data.([vars_names{ii} '_quality_control'])=ncread([url_data glider_type names{i} '/' nc_filename{1}],([vars_names{ii} '_quality_control']));            
        
                %% QC1 - Only 'good data' from the .nc QC-flag
                disp('Loading filtering only "good data"...')
                ind_qc=data.([vars_names{ii} '_quality_control'])~=1;
                data.(vars_names{ii})(ind_qc)=NaN;
    
                clear ind_qc      

                %% QC2 - Bathymetry
                data.(vars_names{ii})(~flags_bathy)=NaN;
                
                %% QC3 - Bottom profile 
                % 20m above topography (bottom boundary layer)
                data.(vars_names{ii})(flags_bbl)=NaN;

                %% QC4 and QC5 - Outliers and quenching on chlorophyll
                % Calculate mean and std over a moving window size of 1000 points
                % and remove points over 2*std
                
                if strcmp(vars_names{ii},'CPHL')  % Check if the variable is chlorophyll
    
                    disp('Performing QCs in chlorophyll...')

                    %% QC4 - Outliers on chlorophyll
                    % Calculate mean and std over a moving window size of 1000 points
                    % and remove points over 2*std
                
                    x=log10(data.(vars_names{ii}));
    
                    mov_mean=movmean(x,window,'omitnan');
                    mov_std=movstd(x,window,'omitnan');
    
                    % Create a mask for outliers (above local mean + 2*std)
                    outlier_mask=find(x>(mov_mean+2*mov_std));
    
                    % Apply mask to set outliers to NaN
                    x_QC=x;
                    x_QC(outlier_mask)=NaN;
    
                    % Optionally smooth with median filter
                    x_smooth=medfilt1(x_QC,5,'omitnan','truncate');
    
                    % Convert back to linear space
                    data.(vars_names{ii})=10.^x_smooth;

                    %% QC5 - Quenching
                    data.(vars_names{ii})(flags_daytime)=NaN;
                    
                end
    
                %% Binning data to fixed depths
                disp('Binning data to fixed depths...')
                cast=casts(cast,data.(vars_names{ii}),vars_names{ii},dims.DEPTH, ...
                           'down',ind_pu,ind_pd,max_depth,depth_res);
                cast=casts(cast,data.(vars_names{ii}),vars_names{ii},dims.DEPTH, ...
                           'up',ind_pu,ind_pd,max_depth,depth_res);
    
                %% Average down-up casts
                disp('Averaging down-up casts for variable...')
                var_dummy=cat(3,cast.down.(vars_names{ii}),cast.up.(vars_names{ii}));
                cast.new.(vars_names{ii})=mean(var_dummy,3,'omitnan');
    
                %% Clean empty coordinates
                cast.new.(vars_names{ii})(is_nan)=[];

                %% Binning data to fixed distances
                disp('Binning data to fixed distances...')
                DATA.(vars_names{ii})=[];
                for iii=1:(max_depth+1)
                    [~,var_dummy]=binning(cast.new.DISTANCE,cast.new.(vars_names{ii})(iii,:), ...
                                          max_dist,dist_res);
                    DATA.(vars_names{ii})=cat(1,DATA.(vars_names{ii}),var_dummy');
                end
    
                clear var_dummy                
    
            catch me
                fprintf('%s\n',me.message);
                fprintf('File: %s\n\n',names{i});
    
                DATA.(vars_names{ii})=NaN(size(DATA.(vars_names{1})));
    
                infos_missing.(names{i}).(vars_names{ii})=1;
    
            end
    
        end

        %% Cleaning...
        % Following temperature variable as reference
        cleaning=mean(DATA.TEMP,1,'omitnan');
        is_nan=isnan(cleaning);
        
        fields=fieldnames(DATA);
        fields(ismember(fields,'DEPTH'))=[];
        for ii=1:length(fields)
            f=fields{ii};
            n_dim=find(size(DATA.(f))==length(cleaning));
            if n_dim==1
                DATA.(f)(is_nan)=[];
            elseif n_dim==2
                DATA.(f)(:,is_nan)=[];
            end
        end

        clear cleaning is_nan f n_dim ii fields

        %% Save DATA in .nc file
        disp(' ')
            
        % Defining variable dimensions and attributes
        NPROF=1:length(DATA.DISTANCE);
    
        if ~isempty(NPROF)
            disp('Saving data in NetCDF file...')
            ncid=netcdf.create([path_nc names{i} '_QCs.nc'],'NC_WRITE');
    
            depthID=netcdf.defDim(ncid,'DEPTH',length(DATA.DEPTH));
            vdepthID=netcdf.defVar(ncid,'DEPTH','float',depthID);
            netcdf.putAtt(ncid,vdepthID,'long_name','Depth');
            netcdf.putAtt(ncid,vdepthID,'units','meters');
        
            nprofID=netcdf.defDim(ncid,'NPROF',length(NPROF));
            vnprofID=netcdf.defVar(ncid,'NPROF','float',nprofID);
            netcdf.putAtt(ncid,vnprofID,'long_name','Profile number');
            netcdf.putAtt(ncid,vnprofID,'units','none');
    
            % Defining variables and attributes
            lonID=netcdf.defVar(ncid,'LONGITUDE','float',nprofID);
            netcdf.putAtt(ncid,lonID,'long_name','Longitude');
            netcdf.putAtt(ncid,lonID,'units','degrees_east');
        
            latID=netcdf.defVar(ncid,'LATITUDE','float',nprofID);
            netcdf.putAtt(ncid,latID,'long_name','Latitude');
            netcdf.putAtt(ncid,latID,'units','degrees_north');
        
            distID=netcdf.defVar(ncid,'DISTANCE','float',nprofID);
            netcdf.putAtt(ncid,distID,'long_name','Distance along-track');
            netcdf.putAtt(ncid,distID,'units','meters');
        
            timeID=netcdf.defVar(ncid,'TIME','float',nprofID);
            netcdf.putAtt(ncid,timeID,'long_name','Time');
            netcdf.putAtt(ncid,timeID,'units','days since 1950-01-01 00:00:00 UTC');
        
            flds=fieldnames(DATA);
            flds_ind=ismember(vars_names,flds);
            flds=vars_names(flds_ind);
            flds_attrb=vars_attrb(flds_ind);
            flds_units=vars_units(flds_ind);
        
            varID=[];
            for ii=1:length(flds)
                varID(ii)=netcdf.defVar(ncid,flds{ii},'float',[depthID nprofID]);
                netcdf.putAtt(ncid,varID(ii),'long_name',flds_attrb{ii});
                netcdf.putAtt(ncid,varID(ii),'units',flds_units{ii});
            end
        
            % End define mode
            netcdf.endDef(ncid)
            
            % Input data
            netcdf.putVar(ncid,vdepthID,DATA.DEPTH);
            netcdf.putVar(ncid,vnprofID,NPROF);
            netcdf.putVar(ncid,lonID,DATA.LONGITUDE);
            netcdf.putVar(ncid,latID,DATA.LATITUDE);
            netcdf.putVar(ncid,distID,DATA.DISTANCE);
            netcdf.putVar(ncid,timeID,DATA.TIME);
            
            for ii=1:length(flds)
                netcdf.putVar(ncid,varID(ii),DATA.(flds{ii}));
            end
        
            % Close the file
            netcdf.close(ncid)
        else
            disp('Empty file after QCs... Skipping!')
        end
    else
        disp('File already exists... Skipping!')
    end

    %% Progress bar for loop of every glider name
    waitbar(i/length(names),pb,sprintf('Progress: %d %%',floor(i/length(names)*100)));
    pause(.01)
    disp(' ')

end

close(pb)

%% Clear workspace
clear i ii iii pb glider_type me nc_filename dist ind_pd ind_pu pks_down pks_up

%% Functions
function [new_reference,new_var]=binning(reference,var,max_ref,ref_res)
    
    new_var=[];
    dummy_depth=(-ref_res/2):ref_res:(max_ref+ref_res/2);
    
    for i=1:(length(dummy_depth)-1)
        ind=reference>=dummy_depth(i) & reference<dummy_depth(i+1);
        new_var=cat(1,new_var,mean(var(ind),'all','omitnan'));
    end

    new_reference=0:ref_res:max_ref;
    new_reference=new_reference';

end

function cast=casts(cast,data_var,var_name,dims_depth,flag_cast_type, ...
                    index_up_peaks,index_down_peaks,max_depth,depth_res)

    % Initialize variables to concatenate
    cast.(flag_cast_type).DEPTH_ref=[];
    cast.(flag_cast_type).(var_name)=[];
    
    % Loop for each cast
    for i=1:length(index_down_peaks)
        
        % Get indexes for the cast
        switch flag_cast_type
            case 'down'
                indexes=index_up_peaks(i):index_down_peaks(i);
            case 'up'
                indexes=index_down_peaks(i):index_up_peaks(i+1);
        end

        % Select depth and variables for the cast
        depth=dims_depth(indexes);
        var=data_var(indexes);
        
        % Binning data to fixed depths
        if i==1
            [new_depth,new_var]=binning(depth,var,max_depth,depth_res);
        else
            [~,new_var]=binning(depth,var,max_depth,depth_res);
        end

        % Concatenate casts
        cast.(flag_cast_type).DEPTH_ref=cat(2,cast.(flag_cast_type).DEPTH_ref,new_depth);
        cast.(flag_cast_type).(var_name)=cat(2,cast.(flag_cast_type).(var_name),new_var);

    end  

end

