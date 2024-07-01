

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Vincent Guigues.
%Creation date: 15/7/2020.
%Last updated: 21/7/2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: given a list of calls, allocates ambulances to calls using
%Heuristic 3 of the paper "".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calls: a list of structures with struture calls(i) storing the data of call i:

%calls(i).time: time of the call

%calls(i).lat: latitude of the call

%calls(i).long: longitude of the call

%calls(i).ih: index of the hospital where the patient of call i is sent

%calls(i).priority: priority of the call in {0,1,2,3,...}. The lower the
%integer, the more the call is prioritary.

%Calls of priority 0 are of the highest priority and the penalized waiting 
%time is survival_function(t) for t time units until the equipment arrives 
%on the scene of the call and penalized_time(t,1) after the equipment has
%arrived and until arriving to hospital.
%For calls of priority i>=1, the penalized waiting time is penalized_time(t,i)

%calls(i).time_on_scene: time the ambulance stays on the scene of the call

%calls(i).time_equipment: time spent by the moto or drone on the scene if
%call(i).type=0

%calls(i).hosp_needed=1 if the patient of the call needs to go to hospital
%and 0 otherwise.

%calls(i).time_at_hospital: time the ambulance waits at hospital before
%being available for another dispatch.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ambulances is a list of structures with ambulances(i) a structure for the
%data of ambulance i with the following fileds:

%ambulances(i).base.lat is the latitude of the last base ambulance i went to or is currently
%forecast to go.

%ambulances(i).base.long is the longitude of the last base ambulance i went to or is currently
%forecast to go.

%ambulances(i).type is the type of ambulance i. Ambulances of type i can
%attend calls of priorities >=i.

%ambulances(i).free_destination.lat is the latitude of the last hospital
%ambulance i went to or is currently forecast to go.

%ambulances(i).free_destination.long is the longitude of the last hospital
%ambulance i went to or is currently forecast to go.

%ambulance(i).arrival_time_at_f_last_trip: time ambulance i was freed from
%the last hospital visited or time it will be freed from an hospital if it is
%currently going to an hospital. If the call is not sent to hospital this
%is the time the ambulance becomes available after serving the call on the scene

%ambulance(i).arrival_time_at_b_last_trip: arrival time of ambulance i
%at the last base visited or time it will arrive at the next base if it is
%in service

%ambulances(i).arrival_time_at_c_last_trip is the arrival time at the last
%visited call scene or the next visited call scene if the ambulance is going
%to a call scene.

%ambulances(i).call.lat is the latitude of the last call scene ambulance i went to or is currently
%forecast to go.

%ambulances(i).call.long is the longitude of the last call scene ambulance i went to or is currently
%forecast to go.

%ambulances(i).speed is the speed of the ambulance.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_ambulances: number of ambulances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_motos: number of motos

%motos(i).base.lat is the latitude of the last base moto i went to or is currently
%forecast to go.

%motos(i).base.long is the longitude of the last base moto i went to or is currently
%forecast to go.

%motos(i).free_destination.lat is the latitude of the last scene call
%moto i went to or is currently forecast to go.

%motos(i).free_destination.long is the longitude of the last scene call
%moto i went to or is currently forecast to go.

%motos(i).arrival_time_at_f_last_trip: time moto i was freed from
%the last call scene visited or time it will be freed from a call scene if it is
%currently going to a call scene. 

%motos(i).arrival_time_at_b_last_trip: arrival time of moto i
%at the last base visited or time it will arrive at the next base if it is
%in service

%motos(i).number_equipment: quantity of equipments a moto can transport

%motos(i).speed is the speed of the moto.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_drones: number of drones

%drones(i).base.lat is the latitude of the last base moto i went to or is currently
%forecast to go.

%drones(i).base.long is the longitude of the last base moto i went to or is currently
%forecast to go.

%drones(i).free_destination.lat is the latitude of the last scene call
%moto i went to or is currently forecast to go.

%drones(i).free_destination.long is the longitude of the last scene call
%moto i went to or is currently forecast to go.

%drones(i).arrival_time_at_f_last_trip: time moto i was freed from
%the last call scene visited or time it will be freed from a call scene if it is
%currently going to a call scene. 

%drones(i).arrival_time_at_b_last_trip: arrival time of moto i
%at the last base visited or time it will arrive at the next base if it is
%in service

%drones(i).number_equipment: quantity of equipments a drone can transport

%drones(i).speed is the speed of the moto.

%drones(i).altitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hospitals: a list of structures giving the locations of hospitals
%hospitals(i).lat is the the latitude of hospital i
%hospitals(i).long is the the longitude of hospital i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_hospitals is the number of hospitals.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_calls: number of calls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The remaining variables are used to store the trajectories of
%ambulances which are further plotted in discretized time.

%ambulances_times,ambulances_trips,trip_type: cell arrays of size
%(1,nb_ambulances): trip i of ambulance j starts at time
%ambulances_times{1,j}(i) from the location with latitude
%ambulances_trips{1,j}(1,i) and longitude ambulances_trips{1,j}(2,i)
%and ends at time ambulances_times{1,j}(i+1) at location with latitude
%ambulances_trips{1,j}(1,i+1) and longitude ambulances_trips{1,j}(2,i+1).
%The type of this trip is given by integer trip_type{1,j}(i)
%with the following values for the trip type: 1: at base, 2: going to
%the call, 3: stays at the scene of the call, 4: going to hospital,
%5: stays at hospital waiting to be freed, 6: going from hospital to base.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bases(i).lat is the latitude of base i
%bases(i).long is the longitude of base i

%nb_base: number of bases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%waiting_to_equipment(i)=0 if call is does not require equipment and otherwise it is the 
%time needed to bring equipment to the patient on the call scene. 

%waiting_on_scene(i): waiting time between the time of call i and the time the ambulance
%arrives on the scene of call i

%waiting_to_hospital(i): waiting time between the instant the ambulance arrives on the scene
%of call i and the instant the ambulance arrives at hospital to bring the
%patient of call i

%which_ambulance(i): ambulance which took care of call i

%waiting_on_scene_penalized(i)=penalized(calls(i).priority)*waiting_on_scene(i)
%[penalized waiting time between the instant of the call and the instant the ambulance
%arrives on the scene of the call]

%waiting_to_hospital_penalized(i)=penalized(calls(i).priority)*waiting_to_hospital(i)
%[penalized waiting time between the instant the ambulance arrives on the scene of the call
%i and the instant the ambulance arrives at the hospital to bring the patient of call i]

%ambulances_times,ambulances_trips,trip_type: the input variables updated
%with the trajectories of the ambulances to attend the calls.

%calls_end(i) is the instant the patient of call i arrives at hospital.

function [waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_minmax_priorities_closest_available_basent(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalized_time,survival_function,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases,RadiusEarth)

calls_end=zeros(1,nb_calls);
waiting_to_equipment==zeros(1,nb_calls);
waiting_on_scene=zeros(1,nb_calls);
waiting_to_hospital=zeros(1,nb_calls);
waiting_on_scene_penalized=zeros(1,nb_calls);
waiting_to_hospital_penalized=zeros(1,nb_calls);
which_ambulance=zeros(1,nb_calls);

queue_size=0;

%queue contains the queue of calls. We only need to store in
%queue the indexes of the calls (the data of the calls is stored
%in calls).
%queue(i) is the index of i-th call in the queue
%calls are ordered by priority with queue(1) the most prioritary call
queue=[];

%event_call=1 if the event treated is call and 0 otherwise (in which case
%the event is "an ambulance is freed")
event_call=1;

%Current time, updated at time instants corresponding to our 2 types of
%events: 1) calls and 2) ambulance is freed.

current_time=calls(1).time;

%index of the last call treated
index_call=1;

%Number of calls already attended
calls_attended=0;

queue=[];

while (calls_attended<nb_calls)
    if (event_call==1)
        %we add the call index to the queue
        queue=[queue;index_call];
        queue_size=queue_size+1;
    end
    
    %While a call in the queue has not been treated:
    %We compute for each call the quickest ambulance to arrive on the scene
    %We consider the call with the largest waiting time (from time of the
    %call until the time the ambulance arrives on the scene) for which
    %the corresponding ambulance is available
    %We allocate that ambulance to the call.
    %We update arrival times considering this allocation.
    %We consider the next call with the largest waiting time for which the
    %corresponding ambulance is available
    %We allocate that ambulance to the call, and so on....
    
    %mintimes(k) is the minimal waiting time for call with index queue(k)
    mintimes=inf*ones(1,length(queue));
    minpentimes=inf*ones(1,length(queue));
    travelTimes=zeros(length(queue),nb_ambulances);
    travelTimesMotos=zeros(length(queue),nb_motos);
    travelTimesDrones=zeros(length(queue),nb_drones);
    %index_amb(k) is the index of t he ambulance which is the quicked to
    %reach call with index queue(k)
    index_ambs=zeros(1,length(queue));
    index_motos=zeros(1,length(queue));
    index_drones=zeros(1,length(queue));
    amb_type=zeros(1,length(queue));
    mintimesAmb=inf*ones(1,length(queue));
    mintimesDrone=inf*ones(1,length(queue));
    mintimesMotos=inf*ones(1,length(queue));

    for k=1:length(queue)
        if (calls(queue(k)).priority==0)
           for j=1:nb_ambulances
                if (ambulances(j).type==1)
                    if (ambulances(j).arrival_time_at_b_last_trip<=current_time)
                        travelTimes(k,j)=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                    elseif (ambulances(j).arrival_time_at_f_last_trip<=current_time)
                        [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,current_time,ambulances(j).speed);
                        travelTimes(k,j)=travel_time(curr_lat,curr_long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                    end
                    if (ambulances(j).arrival_time_at_f_last_trip<=current_time)
                        if ((travelTimes(k,j)<mintimes(k))||((travelTimes(k,j)==mintimes(k))&&(amb_type(k)<ambulances(j).type)))
                            amb_type(k)=ambulances(j).type;
                            index_ambs(k)=j;
                            mintimesAmb(k)=travelTimes(k,j);
                        end
                    end
                end
            end

            for j=1:nb_motos
                    if (motos(j).arrival_time_at_b_last_trip<=current_time)
                        travelTimesMotos(k,j)=travel_time(motos(j).base.lat,motos(j).base.long,calls(queue(k)).lat,calls(queue(k)).long,motos(j).speed);
                    elseif (motos(j).arrival_time_at_f_last_trip<=current_time)
                        [curr_lat,curr_long]=position_between_origin_destination(motos(j).free_destination.lat,
                        motos(j).free_destination.long,motos(j).base.lat,motos(j).base.long,motos(j).arrival_time_at_f_last_trip,current_time,motos(j).speed);
                        travelTimesMotos(k,j)=travel_time(curr_lat,curr_long,calls(queue(k)).lat,calls(queue(k)).long,motos(j).speed);
                    end
                    if (motos(j).arrival_time_at_f_last_trip<=current_time)
                        if ((travelTimesMotos(k,j)<mintimesMotos(k))||((travelTimesMotos(k,j)==mintimesMotos(k))))
                            index_motos(k)=j;
                            mintimesMotos(k)=travelTimes(k,j);
                        end
                    end
            end

            for j=1:nb_drones
                    if (drones(j).arrival_time_at_b_last_trip<=current_time)
                        travelTimesDrones(k,j)=travel_time_geodesic(drones(j).base.lat,drones(j).base.long,calls(queue(k)).lat,calls(queue(k)).long,drones(j).speed,'spheric',RadiusEarth+drones(j).altitude);
                    elseif (drones(j).arrival_time_at_f_last_trip<=current_time)
                        [curr_lat,curr_long]=position_between_origin_destination_geodesic(drones(j).free_destination.lat,drones(j).free_destination.long,drones(j).base.lat,drones(j).base.long,drones(j).arrival_time_at_f_last_trip,current_time,drones(j).speed);
                        travelTimesDrones(k,j)=travel_time_geodesic(curr_lat,curr_long,calls(queue(k)).lat,calls(queue(k)).long,motos(j).speed,'spheric',RadiusEarth+drones(j).altitude);
                    end
                    if (drones(j).arrival_time_at_f_last_trip<=current_time)
                        if ((travelTimesDrones(k,j)<mintimesDrones(k))||((travelTimesDrones(k,j)==mintimesDrones(k))))
                            index_drones(k)=j;
                            mintimesDrones(k)=travelTimesDrones(k,j);
                        end
                    end
            end

            if ((mintimesMotos(k)<mintimesAmb(k))||(mintimesDrones(k)<mintimesAmb(k)) 
                 mintimes(k)=mintimesAmb(k);
                 if (mintimesMotos(k)<mintimesDrones(k))
                     minpentimes(k)=survival_function(mintimesMotos(k))+penalized_time(mintimesAmb(k)-mintimesMotos(k),1);
                     index_drones(k)=0;            
                 else
                     minpentimes(k)=survival_function(mintimesDrones(k))+penalized_time(mintimesAmb(k)-mintimesDrones,1);
                     index_motos(k)=0;
                 end 
            else
               mintimes(k)=mintimesAmb(k);
               minpentimes(k)=survival_function(mintimes(k));
               index_drones(k)=0;
               index_motos(k)=0;                   
            end
        else
            for j=1:nb_ambulances
                if (ambulances(j).type<=calls(queue(k)).priority)
                    if (ambulances(j).arrival_time_at_b_last_trip<=current_time)
                        travelTimes(k,j)=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                    elseif (ambulances(j).arrival_time_at_f_last_trip<=current_time)
                        %Where is the ambulance?
                        %It started at ambulance(j).arrival_time_at_f_last_trip its
                        %travel from [ambulances(j).free_destination.lat,ambulances(j).free_destination.long]
                        %to its base [ambulances(j).base.lat,ambulances(j).base.long]
                        %Where is it at current_time?
                        [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,current_time,ambulances(j).speed);
                        travelTimes(k,j)=travel_time(curr_lat,curr_long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                    end
                    if (ambulances(j).arrival_time_at_f_last_trip<=current_time)
                        if ((travelTimes(k,j)<mintimes(k))||((travelTimes(k,j)==mintimes(k))&&(amb_type(k)<ambulances(j).type)))
                            amb_type(k)=ambulances(j).type;
                            index_ambs(k)=j;
                            mintimes(k)=travelTimes(k,j);
                        end
                    end
                end
            end
        end
        mintimes(k)=mintimes(k)+current_time-calls(queue(k)).time;
        minpentimes(k)=penalized_time(calls(queue(k)).priority)*mintimes(k);
    end
    
    remainingIndexes=[1:length(queue)];
    
    queue_aux=[];
    laux=0;
    nbTreateadCalls=0;
    totalInQueue=length(queue);
    
    while (nbTreateadCalls<totalInQueue)
        %We treat call with index queue(remainingIndexes(currentCall)) and the best ambulance is with index index_ambs(remainingIndexes(currentCall))
        %[thismintimes,currentCall]=max(mintimes);
        [thispenmintimes,currentCall]=max(minpentimes);
        thismintimes=mintimes(currentCall);
        indexCall=queue(remainingIndexes(currentCall));
        index_amb=index_ambs(remainingIndexes(currentCall));
        if (index_amb>0)
            if (ambulances(index_amb).arrival_time_at_f_last_trip<=current_time)
                if (ambulances(index_amb).arrival_time_at_b_last_trip<=current_time)
                    ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},current_time];
                    ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[ambulances(index_amb).base.lat;ambulances(index_amb).base.long]];
                    trip_type{1,index_amb}=[trip_type{1,index_amb},1];
                else
                    if  (ambulances(index_amb).arrival_time_at_f_last_trip==current_time)
                        len=length(ambulances_times{1,index_amb});
                        ambulances_trips{1,index_amb}=ambulances_trips{1,index_amb}(:,1:len-1);
                        ambulances_times{1,index_amb}=ambulances_times{1,index_amb}(:,1:len-1);
                        trip_type{1,index_amb}=trip_type{1,index_amb}(:,1:len-2);
                    else
                        len=length(ambulances_times{1,index_amb});
                        [curr_lat,curr_long]=position_between_origin_destination(ambulances(index_amb).free_destination.lat,ambulances(index_amb).free_destination.long,ambulances(index_amb).base.lat,ambulances(index_amb).base.long,ambulances(index_amb).arrival_time_at_f_last_trip,current_time,ambulances(index_amb).speed);
                        ambulances_trips{1,index_amb}(:,len)=[curr_lat;curr_long];
                        ambulances_times{1,index_amb}(len)=current_time;
                    end
                end
                
                timeArrivalScene=calls(indexCall).time+thismintimes;
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeArrivalScene];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[calls(indexCall).lat;calls(indexCall).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},2];
                
                timeLeaveScene=timeArrivalScene+calls(indexCall).time_on_scene;
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeLeaveScene];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[calls(indexCall).lat;calls(indexCall).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},3];
                
                calls_attended=calls_attended+1;
                which_ambulance(indexCall)=index_amb;
                waiting_on_scene(indexCall)=thismintimes;
                
                if (calls(indexCall).cleaning_needed)
                    if (calls(indexCall).hosp_needed)
                        t_from_c_to_h=travel_time(calls(indexCall).lat,calls(indexCall).long,hospitals(calls(indexCall).ih).lat,hospitals(calls(indexCall).ih).long,ambulances(index_amb).speed);
                        indexcb=hospitals(calls(indexCall).ih).cbase;
                        clbase=cleaning_bases(indexcb);
                        finalBase=cleaning_bases(indexcb).base;
                        t_from_h_to_cb=travel_time(hospitals(calls(indexCall).ih).lat,hospitals(calls(indexCall).ih).long,clbase.lat,clbase.long,ambulances(index_amb).speed);
                        t_from_cb_to_b=travel_time(clbase.lat,clbase.long,finalBase.lat,finalBase.long,ambulances(index_amb).speed);
                        tLeaveH=timeLeaveScene+calls(indexCall).time_at_hospital+t_from_c_to_h;
                        tArrivalcb=tLeaveH+t_from_h_to_cb;
                        tLeavecb=tArrivalcb+calls(indexCall).cleaning_time;
                        tBase=tLeavecb+t_from_cb_to_b;
                        waiting_to_hospital(indexCall)=calls(indexCall).time_on_scene+t_from_c_to_h;
                        ambulances(index_amb).free_destination.lat=clbase.lat;
                        ambulances(index_amb).free_destination.long=clbase.long;
                        
                        ambulances(index_amb).arrival_time_at_f_last_trip=tLeavecb;
                        ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
                        ambulances(index_amb).base.lat=finalBase.lat;
                        ambulances(index_amb).base.long=finalBase.long;
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(indexCall).time_at_hospital];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(indexCall).ih).lat;hospitals(calls(indexCall).ih).long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},4];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(indexCall).ih).lat;hospitals(calls(indexCall).ih).long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},5];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tArrivalcb];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[clbase.lat;clbase.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},6];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeavecb];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[clbase.lat;clbase.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},7];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[finalBase.lat;finalBase.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                        
                        calls_end(indexCall)=tLeaveH-calls(indexCall).time_at_hospital;
                    else
                        mindist=inf;
                        indmin=0;
                        for indk=1:length(cleaning_bases)
                            [d]=mydistance(calls(indexCall).lat,calls(indexCall).long,cleaning_bases(indk).lat,cleaning_bases(indk).long);
                            if (d<mindist)
                                mindist=d;
                                indmin=indk;
                            end
                        end
                        cb=cleaning_bases(indmin);
                        t_from_c_to_cb=travel_time(calls(indexCall).lat,calls(indexCall).long,cb.lat,cb.long,ambulances(index_amb).speed);
                        finalBase=cleaning_bases(indmin).base;
                        t_from_cb_to_b=travel_time(cb.lat,cb.long,finalBase.lat,finalBase.long,ambulances(index_amb).speed);
                        timeArrivalcb=timeLeaveScene+t_from_c_to_cb;
                        timeLeavecb=timeArrivalcb+calls(indexCall).cleaning_time;
                        timeBase=timeLeavecb+t_from_cb_to_b;
                        waiting_to_hospital(indexCall)=calls(indexCall).time_on_scene;
                        ambulances(index_amb).free_destination.lat=cb.lat;
                        ambulances(index_amb).free_destination.long=cb.long;
                        ambulances(index_amb).arrival_time_at_f_last_trip=timeLeavecb;
                        ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
                        ambulances(index_amb).base.lat=finalBase.lat;
                        ambulances(index_amb).base.long=finalBase.long;
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeArrivalcb];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cb.lat;cb.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},6];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeLeavecb];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cb.lat;cb.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},7];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeBase];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[finalBase;finalBase]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                        
                        calls_end(indexCall)=timeLeaveScene;
                    end
                else
                    if (calls(indexCall).hosp_needed)
                        t_from_c_to_h=travel_time(calls(indexCall).lat,calls(indexCall).long,hospitals(calls(indexCall).ih).lat,hospitals(calls(indexCall).ih).long,ambulances(index_amb).speed);
                        t_from_h_to_b=travel_time(hospitals(calls(indexCall).ih).lat,hospitals(calls(indexCall).ih).long,hospitals(calls(indexCall).ih).base.lat,hospitals(calls(indexCall).ih).base.long,ambulances(index_amb).speed);
                        tLeaveH=timeLeaveScene+calls(indexCall).time_at_hospital+t_from_c_to_h;
                        tBase=tLeaveH+t_from_h_to_b;
                        waiting_to_hospital(indexCall)=calls(indexCall).time_on_scene+t_from_c_to_h;
                        ambulances(index_amb).free_destination.lat=hospitals(calls(indexCall).ih).lat;
                        ambulances(index_amb).free_destination.long=hospitals(calls(indexCall).ih).long;
                        ambulances(index_amb).arrival_time_at_f_last_trip=tLeaveH;
                        ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
                        ambulances(index_amb).base.lat=hospitals(calls(indexCall).ih).base.lat;
                        ambulances(index_amb).base.long=hospitals(calls(indexCall).ih).base.long;
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(indexCall).time_at_hospital];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(indexCall).ih).lat;hospitals(calls(indexCall).ih).long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},4];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(indexCall).ih).lat;hospitals(calls(indexCall).ih).long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},5];
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[ambulances(index_amb).base.lat;ambulances(index_amb).base.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                        
                        calls_end(indexCall)=tLeaveH-calls(indexCall).time_at_hospital;
                    else
                        mindist=inf;
                        indmin=0;
                        for indk=1:nb_base
                            [d]=mydistance(calls(indexCall).lat,calls(indexCall).long,bases(indk).lat,bases(indk).long);
                            if (d<mindist)
                                mindist=d;
                                indmin=indk;
                            end
                        end
                        t_from_c_to_b=travel_time(calls(indexCall).lat,calls(indexCall).long,bases(indmin).lat,bases(indmin).long,ambulances(index_amb).speed);
                        tBase=timeLeaveScene+t_from_c_to_b;
                        waiting_to_hospital(indexCall)=calls(indexCall).time_on_scene;
                        
                        ambulances(index_amb).base.lat=bases(indmin).lat;
                        ambulances(index_amb).base.long=bases(indmin).long;
                        ambulances(index_amb).free_destination.lat=calls(indexCall).lat;
                        ambulances(index_amb).free_destination.long=calls(indexCall).long;
                        ambulances(index_amb).arrival_time_at_f_last_trip=timeLeaveScene;
                        ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
                        
                        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
                        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[ambulances(index_amb).base.lat;ambulances(index_amb).base.long]];
                        trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                        
                        calls_end(indexCall)=calls(indexCall).time+mintimes(1)+calls(indexCall).time_on_scene;
                    end
                end
                waiting_on_scene_penalized(indexCall)=penalized_time(waiting_on_scene(indexCall),calls(indexCall).priority);
                waiting_to_hospital_penalized(indexCall)=penalized_time(waiting_to_hospital(indexCall),calls(indexCall).priority);
                
                %Update travelTimes
                for k=1:(totalInQueue-nbTreateadCalls)
                    %Compute mintimes(k) for k~=currentCall which is the minimal waiting time for
                    %call with index queue(remainingIndex(k))
                    if ((k~=currentCall)&&(index_amb==index_ambs(remainingIndexes(k))))
                        %New time for ambulance index_amb to go to call queue(remainingIndexes(k))
                        time_to_hospital=ambulances(index_amb).arrival_time_at_f_last_trip-current_time;
                        time_from_hosp_to_call=travel_time(ambulances(index_amb).free_destination.lat,ambulances(index_amb).free_destination.long,calls(queue(remainingIndexes(k))).lat,calls(queue(remainingIndexes(k))).long,ambulances(index_amb).speed);
                        travelTimes(remainingIndexes(k),index_amb)=time_to_hospital+time_from_hosp_to_call;
                        [value,index]=min(travelTimes(remainingIndexes(k),:));
                        mintimes(k)=value+current_time-calls(queue(remainingIndexes(k))).time;
                        index_ambs(remainingIndexes(k))=index;
                        minpentimes(k)=penalized_time(mintimes(k),calls(queue(remainingIndexes(k))).priority);
                    end
                end
            else
                queue_aux=[queue_aux;queue(remainingIndexes(currentCall))];
                laux=laux+1;
            end
        else
            queue_aux=[queue_aux;queue(remainingIndexes(currentCall))];
            laux=laux+1;
        end
        if ((nbTreateadCalls+1)<totalInQueue)
            mintimes=[mintimes(1:currentCall-1);mintimes(currentCall+1:totalInQueue-nbTreateadCalls)];
            remainingIndexes=[remainingIndexes(1:currentCall-1);remainingIndexes(currentCall+1:totalInQueue-nbTreateadCalls)];
            minpentimes=[minpentimes(1:currentCall-1);minpentimes(currentCall+1:totalInQueue-nbTreateadCalls)];
        end
        nbTreateadCalls=nbTreateadCalls+1;
    end
    
    queue_size=laux;
    queue=queue_aux;
    
    %Update current_time and index_call
    %future_arrival_time_hospital gives the future time instant where
    %ambulances will arrive at some hospitals
    future_arrival_time_hospital=[];
    for j=1:nb_ambulances
        if (ambulances(j).arrival_time_at_f_last_trip>current_time)
            future_arrival_time_hospital=[future_arrival_time_hospital;ambulances(j).arrival_time_at_f_last_trip];
        end
    end
    
    if (length(future_arrival_time_hospital)>0)
        minh=min(future_arrival_time_hospital);
        if (index_call<nb_calls)
            if (calls(index_call+1).time<=minh)
                event_call=1;
                index_call=index_call+1;
                current_time=calls(index_call).time;
            else
                event_call=0;
                current_time=minh;
            end
        else
            event_call=0;
            current_time=minh;
        end
    else
        if (index_call<nb_calls)
            event_call=1;
            index_call=index_call+1;
            current_time=calls(index_call).time;
        end
    end
end





