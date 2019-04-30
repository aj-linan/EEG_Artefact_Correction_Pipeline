% 02-04-2019 
% Author: Albert Linan
% Instituto de Sistemas e Robotica Lisboa
% This funcion creates the markers for deleting the GA(Gradient Artifact)
% using the Time of Repetition of a slice(TRslice) and a Threshold(TH)

function [new_alldata, new_data, points] = create_markers_v6(alldata, data, TR , TH, figures, s_channels, volumes)

    %% Inputs and outputs
    
        % Inputs

        % alldata ~ ALLEEG struct
        % data ~ EEG struct

        % Time of repetition(ms)
        % TR = 1260; 

        % Threshold
        % This threshold should be selected loocking at the histogram but 200 works
        % TH = 200;

        % Figures
        % 1: Displays the histogram and the rest of the figures
        % 0: No 

        % s_channels 
        % It is a vector that contains the selected channels to calculate the 
        % average channel used for creating the markers

        % volumes
        % 1: Create volume/slice markers
        % 0: Create slices markers

    % Outputs
    
        % new_alldata ~ ALLEEG struct
        % new_data ~ EEG struct
        % points ~ contains all the markers
    
    warning('off','all')

    %% Lets find the points where the gradient grows fast
       
    new_data = data;
    new_alldata = alldata;
    
    % Sampling rate(Hz)
    srate = data.srate; 
        
    % Time of repetition(points)
    TRpoints = TR*srate/1000; 
    
     % Number of slices per volume. Usually 20
    nslicevol = 20;
    
    % Time of repetition for each slice(points)
    TRpointslice = TRpoints/nslicevol;
    
    % If the field 's_channels' is empty, the average of all channels is calculated
    % Else, average is performed with the selected channels
    
    if isempty(s_channels)
        
        disp('All channels selected')
        m_channel = mean(data.data);
    
    else
        
        n_channel = data.nbchan;
        l_channel = 1:n_channel;
        template = ismember(l_channel,s_channels);
        m_channel = (template * [data.data]  /sum(template));
        
    end
    
    % Ploting the histogram
    if figures
        
        figure
        histogram(abs(gradient(m_channel)),20)
        title('Gradient Histogram')
        
    end
    
    % The difference is calculated and the threshold is applied 
    % Y = diff(X)
    % Y = [X(2)-X(1) X(3)-X(2) ... X(m)-X(m-1)]
    
    ids = find(diff(m_channel) > TH);
    
    % Threshold for selecting the good ones
    
    THcont = median(unique(diff(ids))) - 5;
    
    ids_slice = find(diff(ids)> THcont);
    
    ids_slice2 = find(abs(gradient(ids)) > TRpoints);

    % If there are not any marker created with the threshold selected the function ends
    % Else, the new structure with the markers is created
    
    if isempty(ids)
        
        disp('No markers created, try another threshold')
        points = [];
    
    else
    
        % The vector 'points' contains all the markers
        % The vector 'points2' contains the border of the gradient

        points = ids([ids_slice,length(ids)]);
        
        points2 = ids([ids_slice2, length(ids)]);
        
        % The points are shifted to the beggining of each artifact
        
        dist = diff(points);
        dist = dist(1:10);
        mediana = median(dist);
        f_point = find(dist >( mediana + 10));
        
        if points2(1) <= points(1)
        
            if ~isempty(f_point)
                
                id = points2(f_point(1)+1);
                diferencia = points(f_point +1) - id;
                points = points - diferencia.*ones(size(points)); 
                
            else
                
                id = points2(1);
                diferencia = points(f_point) - id;
                points = points - diferencia.*ones(size(points));
                    
            end
        
        end   
       
        % points with diferente TR are stored
        vol = find(diff(points) > 1.1 * TRpointslice & diff(points) < 2 * TRpointslice);
        
        if isempty(vol)
            
            disp('No markers created, try another threshold or check the TR')
            points = [];
        
        else  
            
            vol = [vol(1) - (nslicevol - 1), vol + 1];
            points_vol = points(vol);

            dummy_points = [];
            rest_points = [];

            if ~isempty(points_vol)

                dummy_points = points(points < points_vol(1)); 
                rest_points = points(points > points_vol(1)); 

            end

            points = points(2: length(points));
            points = points(points >= points_vol(1));

            %% Displays
            % Number of slices or volumes/slice are displayed

            disp1 = [' Have been found ', num2str(length(points)) ,' slices. ' ];
            disp(disp1);

            disp2 = [' Have been found ', num2str(length(points_vol)), ' volumes. '];
            disp(disp2);

            if figures

                figure
                subplot(2,2,[3,4]);
                hold on
                plot((m_channel))

                plot(points,zeros(length(points),1),'r*')    
                plot(points_vol,zeros(length(points_vol),1),'go')

                legend('Channels Average', 'Slices' , 'Volumes')
                xlabel('Points')
                ylabel('Amplitude(uV)')
                title('Full averaged signal with markers')
                hold off

                subplot(2,2,1);
                histogram(diff(points),5);
                title('Slices Histogram')

                subplot(2,2,2);
                histogram(diff(points_vol),5);
                title('Volumes Histogram')

            end       

            % Depepending on the format data.event(1,1) is empty or [], for
            % this code must be []. It will not affect the results neither the
            % operations in 'eeglab'
            data.event(1,1).type = [];

            %% Calculation of the rest of the markers
            % The fields of each event follow the same structure as 'eeglab' files
            % The value of the fields Type and Value are the same that
            % 'Analyzer uses

            % 'Latency'
            % Choose working with volumes/slices or slices

            if volumes 

                event_latency1 = points_vol;    

            else

                event_latency1 = points( points >= points_vol(1) );

            end

            n_events = length(event_latency1);
            event_latency = num2cell([event_latency1, data.event.latency]);

            % 'type'. All the markers has the same type: 'Scan Start' 
            event_type = strsplit(repmat('Scan Start,',1,n_events),',');
            event_type1 = event_type(1:n_events);
            event_type = [event_type1, ' ' , data.event.type];

            % Thes comment depends on the eeglab version
%             % 'value'. All the markers has the same value: 'Scanner'
%             event_value = strsplit(repmat('Scanner,',1,n_events),',');
%             event_value = event_value(1:n_events);
%             event_value = [event_value, data.event.value];

            % 'duration'
            event_duration = num2cell([ones(1, n_events),data.event.duration]);

            % The estructure is created

            field1 = 'type'; value1 = event_type;
            % This comment depends on the eeglab version
%             field2 = 'value'; value2 = event_value;
            field3 = 'latency'; value3 = event_latency;
            field4 = 'duration'; value4 = event_duration;

            % This comment depends on the eeglab version
%             own_event = struct(field1, value1, field2, value2, field3, value3, field4, value4);

            own_event = struct(field1, value1, field3, value3, field4, value4);

            [~,index] = sortrows({own_event.latency}.'); 
            own_event = own_event(index);

            % The field 'urevent' is added to the structure
            event_id = num2cell(1:length(own_event));
            [own_event(:).urevent] = deal(event_id{:});

            % Acording with the eeglab files structures, It is needed to create the
            % struct urevent
            own_urevent = rmfield(own_event,'urevent');

            %% The structures are updated

            points = event_latency1;

            new_data.event = own_event;
            new_data.urevent = own_urevent;

            n_data = length(alldata);
            new_alldata(n_data + 1) = new_data;
            
        end
        
    end
end