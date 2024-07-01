 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Vincent Guigues.
%Creation date: 14/7/2020.
%Last updated: 21/7/2020.
%Validated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: given a list of calls, allocates ambulances to calls using
%Heuristic Best myopic of the paper "New Heuristics for the Operation of an Ambulance 
%Fleet under Uncertainty" by V. Guigues, A. Kleywegt, V.H. Nascimento.
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

%calls(i).priority: priority of the call in {1,2,3,...}. The lower the
%integer, the higher the priority.

%calls(i).time_on_scene: time the ambulance stays on the scene of the call

%calls(i).hosp_needed=1 if the patient of the call needs to go to hospital
%and 0 otherwise.

%calls(i).cleaning_needed=1 if the ambulance in charge of call i needs 
%cleaning after service and 0 otherwise.

%calls(i).time_at_hospital: time the ambulance waits at hospital before
%being available for another dispatch.

%calls(i).cleaning_time: time required to clean the ambulance attending
%call i.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ambulances is a list of structures with ambulances(i) a structure for the
%data of ambulance i with the following fields: 

%ambulances(i).base.lat is the latitude of the last base ambulance i went to or is currently 
%forecast to go.

%ambulances(i).base.long is the longitude of the last base ambulance i went to or is currently 
%forecast to go.

%ambulances(i).type is the type of ambulance i. Ambulances of type i can
%attend calls of priorities >=i. 

%ambulances(i).free_destination.lat is the latitude of the last location
%ambulance i was freed or is currently forecast to be freed.

%ambulances(i).free_destination.long is the longitude of the last location
%ambulance i was freed or is currently forecast to be freed.

%ambulance(i).arrival_time_at_f_last_trip: time ambulance i was freed from
%the last location it became available or time it will become available if it is
%currently in service.

%ambulance(i).arrival_time_at_b_last_trip: arrival time of ambulance i
%at the last base visited or time it will arrive at the next base if it is
%in service

%ambulances(i).speed is the speed of the ambulance.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_ambulances: number of ambulances 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hospitals: a list of structures giving the locations of hospitals
%hospitals(i).lat is the the latitude of hospital i
%hospitals(i).long is the the longitude of hospital i
%hospitals(i).cbase.lat is the latitude of the closest cleaning base to hospital i
%hospitals(i).cbase.long is the longitude of the closest cleaning base to hospital i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_hospitals is the number of hospitals.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_calls: number of calls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%penalization: Used to compute penalized waiting times.
%For calls of priority i, the penalized waiting time is  
%(Real waiting time)*penalization(i).

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
%5: stays at hospital waiting to be freed, 6: going to a cleaning base
%7: the ambulance is being cleaned at a base,
%8: the ambulance is going to a base (not for a cleaning task).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bases(i).lat is the latitude of base i
%bases(i).long is the longitude of base i

%nb_base: number of bases

%list_cleaning_bases(i) is the index of ith cleaning base.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_forward_closest_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases)

calls_end=zeros(1,nb_calls);
waiting_on_scene=zeros(1,nb_calls);
waiting_to_hospital=zeros(1,nb_calls);
waiting_on_scene_penalized=zeros(1,nb_calls);
waiting_to_hospital_penalized=zeros(1,nb_calls);
which_ambulance=zeros(1,nb_calls);

for i=1:nb_calls
    mintime=inf;
    index_amb=0;
    %Computation of the quickest response time for call i given the previous
    %allocation of ambulances to calls 1,2,...,i-1.  
    
    for j=1:nb_ambulances
        %First check if ambulance j can take calls of priority calls(i).priority)
        if (ambulances(j).type<=calls(i).priority)
            %Computation of the time for ambulance j to go to the scene of call i, i.e., [calls(i).lat,calls(i).long]
            %This time is stored in time_this_ambulance
            %We first need to know when ambulance j will be free
            %Is it already available when call i arrives? This is the case iff
            %(a) it is at its base when this call arrives or
            %(b) it is between an hospital and its base, going to this base, when the call arrives
            %(a) is equivalent to
            %ambulance(j).arrival_time_at_b_last_trip<=calls(i).time
            %(b) is equivalent to
            %ambulance(j).arrival_time_at_b_last_trip>calls(i).time
            %and ambulance(j).arrival_time_at_f_last_trip<=calls(i).time
            %In the first case it will go from its base to the call
            %In the second case, it will go from its current position (on its way from the hospital
            %to the base) to the call
            if (ambulances(j).arrival_time_at_b_last_trip<=calls(i).time)
                time_this_ambulance=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(i).lat,calls(i).long,ambulances(j).speed);
            elseif (ambulances(j).arrival_time_at_f_last_trip<=calls(i).time)
                %where is the ambulance?
                %It started at ambulance(j).arrival_time_at_f_last_trip its
                %ride from [ambulances(j).free_destination.lat,ambulances(j).free_destination.long]
                %to its base [ambulances(j).base.lat,ambulances(j).base.long]
                %Computation of its location at calls(i).time      
                [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,calls(i).time,ambulances(j).speed);
                time_this_ambulance=travel_time(curr_lat,curr_long,calls(i).lat,calls(i).long,ambulances(j).speed);
            else
                %The ambulance is not free if and only if
                %ambulance(j).arrival_time_at_f_last_trip>calls(i).time
                %In this case to go to [calls(i).lat,calls(i).long], the ambulance will
                %first go to the location it is currently headed to,
                %waiting to be freed, and go to the call.
                %time_free is the time needed before becoming available
                time_free=ambulances(j).arrival_time_at_f_last_trip-calls(i).time;
                time_from_free_to_call=travel_time(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,calls(i).lat,calls(i).long,ambulances(j).speed);
                time_this_ambulance=time_free+time_from_free_to_call;
            end
            if ((time_this_ambulance<mintime)||((time_this_ambulance==mintime)&&(ambulances(j).type>ambulances(index_amb).type)))
                index_amb=j;
                mintime=time_this_ambulance;
            end
        end
    end
              
    %If index_amb=0 this means that the call cannot be attended because no
    %ambulance can attend calls of the corresponding priority.
    %Such situation should not happen in practise.
    if (index_amb>0)        
        if (ambulances(index_amb).arrival_time_at_b_last_trip<calls(i).time)
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},calls(i).time];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[ambulances(index_amb).base.lat;ambulances(index_amb).base.long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},1];
        elseif (ambulances(index_amb).arrival_time_at_f_last_trip<=calls(i).time)
            if  (ambulances(index_amb).arrival_time_at_f_last_trip==calls(i).time)
                 len=length(ambulances_times{1,index_amb});
                 ambulances_trips{1,index_amb}=ambulances_trips{1,index_amb}(:,1:len-1);
                 ambulances_times{1,index_amb}=ambulances_times{1,index_amb}(:,1:len-1);
                 trip_type{1,index_amb}=trip_type{1,index_amb}(:,1:len-2);
            else
                len=length(ambulances_times{1,index_amb});
                [curr_lat,curr_long]=position_between_origin_destination(ambulances(index_amb).free_destination.lat,ambulances(index_amb).free_destination.long,ambulances(index_amb).base.lat,ambulances(index_amb).base.long,ambulances(index_amb).arrival_time_at_f_last_trip,calls(i).time,ambulances(index_amb).speed);
                ambulances_trips{1,index_amb}(:,len)=[curr_lat;curr_long];
                ambulances_times{1,index_amb}(len)=calls(i).time;
            end
        else
            len=length(ambulances_times{1,index_amb});
            ambulances_trips{1,index_amb}=ambulances_trips{1,index_amb}(:,1:len-1);
            ambulances_times{1,index_amb}=ambulances_times{1,index_amb}(:,1:len-1);
            trip_type{1,index_amb}=trip_type{1,index_amb}(:,1:len-2);
        end
        
        timeArrivalScene=calls(i).time+mintime;
        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeArrivalScene];
        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[calls(i).lat;calls(i).long]];
        trip_type{1,index_amb}=[trip_type{1,index_amb},2];
        
        timeLeaveScene=timeArrivalScene+calls(i).time_on_scene;
        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeLeaveScene];
        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[calls(i).lat;calls(i).long]];
        trip_type{1,index_amb}=[trip_type{1,index_amb},3];
        
        which_ambulance(i)=index_amb;
        waiting_on_scene(i)=mintime;
        
        if (calls(i).cleaning_needed)
            if (calls(i).hosp_needed)
                t_from_c_to_h=travel_time(calls(i).lat,calls(i).long,hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,ambulances(index_amb).speed);
                indexcb=hospitals(calls(i).ih).cbase;
                cbase=cleaning_bases(indexcb);
                t_from_h_to_cb=travel_time(hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,cbase.lat,cbase.long,ambulances(index_amb).speed);
                t_from_cb_to_b=travel_time(cbase.lat,cbase.long,cleaning_bases(indexcb).base.lat,cleaning_bases(indexcb).base.long,ambulances(index_amb).speed);
                tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
                tArrivalcb=tLeaveH+t_from_h_to_cb;
                tLeavecb=tArrivalcb+calls(i).cleaning_time;
                tBase=tLeavecb+t_from_cb_to_b;
                
                calls_end(i)=timeLeaveScene+t_from_c_to_h;
                waiting_to_hospital(i)=calls(i).time_on_scene+t_from_c_to_h;
                
                ambulances(index_amb).free_destination.lat=cbase.lat;
                ambulances(index_amb).free_destination.long=cbase.long;
                ambulances(index_amb).arrival_time_at_f_last_trip=tLeavecb;
                ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
                ambulances(index_amb).base.lat=cleaning_bases(indexcb).base.lat;
                ambulances(index_amb).base.long=cleaning_bases(indexcb).base.long;
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(i).time_at_hospital];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},4];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
                 ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},5];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tArrivalcb];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cbase.lat;cbase.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},6];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeavecb];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cbase.lat;cbase.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},7];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cleaning_bases(indexcb).base.lat;cleaning_bases(indexcb).base.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            else                
                mindist=inf;
                indmin=0;
                for indk=1:length(cleaning_bases)
                    [d]=mydistance(calls(i).lat,calls(i).long,cleaning_bases(indk).lat,cleaning_bases(indk).long);
                    if (d<mindist)
                        mindist=d;
                        indmin=indk;
                    end
                end
                cb=cleaning_bases(indmin);
                thisBase=cleaning_bases(indmin).base;
                t_from_c_to_cb=travel_time(calls(i).lat,calls(i).long,cb.lat,cb.long,ambulances(index_amb).speed);
                t_from_cb_to_b=travel_time(cb.lat,cb.long,thisBase.lat,thisBase.long,ambulances(index_amb).speed);
                timeArrivalcb=timeLeaveScene+t_from_c_to_cb;
                timeLeavecb=timeArrivalcb+calls(i).cleaning_time;    
                timeBase=timeLeavecb+t_from_cb_to_b;
                waiting_to_hospital(i)=calls(i).time_on_scene;
                ambulances(index_amb).free_destination.lat=cb.lat;
                ambulances(index_amb).free_destination.long=cb.long;
                ambulances(index_amb).arrival_time_at_f_last_trip=timeLeavecb;
                ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
                ambulances(index_amb).base.lat=thisBase.lat;
                ambulances(index_amb).base.long=thisBase.long;
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeArrivalcb];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cb.lat;cb.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},6];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeLeavecb];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[cb.lat;cb.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},7];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeBase];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[thisBase.lat;thisBase.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                
                calls_end(i)=timeLeaveScene;
            end
        else            
            if (calls(i).hosp_needed)
                t_from_c_to_h=travel_time(calls(i).lat,calls(i).long,hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,ambulances(index_amb).speed);
                t_from_h_to_b=travel_time(hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,hospitals(calls(i).ih).base.lat,hospitals(calls(i).ih).base.long,ambulances(index_amb).speed);
                tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
                tBase=tLeaveH+t_from_h_to_b;
                waiting_to_hospital(i)=calls(i).time_on_scene+t_from_c_to_h;
                ambulances(index_amb).free_destination.lat=hospitals(calls(i).ih).lat;
                ambulances(index_amb).free_destination.long=hospitals(calls(i).ih).long;
                ambulances(index_amb).arrival_time_at_f_last_trip=tLeaveH;
                ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
                ambulances(index_amb).base.lat=hospitals(calls(i).ih).base.lat;
                ambulances(index_amb).base.long=hospitals(calls(i).ih).base.long;
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(i).time_at_hospital];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},4];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},5];
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).base.lat;hospitals(calls(i).ih).base.long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                calls_end(i)=timeLeaveScene+t_from_c_to_h;
            else                
                mindist=inf;
                indmin=0;
                for indk=1:nb_base
                    [d]=mydistance(calls(i).lat,calls(i).long,bases(indk).lat,bases(indk).long);
                    if (d<mindist)
                        mindist=d;
                        indmin=indk;
                    end
                end
                t_from_c_to_b=travel_time(calls(i).lat,calls(i).long,bases(indmin).lat,bases(indmin).long,ambulances(index_amb).speed);
                timeBase=timeLeaveScene+t_from_c_to_b;
                waiting_to_hospital(i)=calls(i).time_on_scene;
                ambulances(index_amb).free_destination.lat=calls(i).lat;
                ambulances(index_amb).free_destination.long=calls(i).long;
                ambulances(index_amb).arrival_time_at_f_last_trip=timeLeaveScene;
                ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
                ambulances(index_amb).base.lat=bases(indmin).lat;
                ambulances(index_amb).base.long=bases(indmin).long;
                
                ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeBase];
                ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[bases(indmin).lat;bases(indmin).long]];
                trip_type{1,index_amb}=[trip_type{1,index_amb},8];
                
                calls_end(i)=timeLeaveScene;
            end
        end
        
        waiting_on_scene_penalized(i)=penalization(calls(i).priority)*waiting_on_scene(i);
        waiting_to_hospital_penalized(i)=penalization(calls(i).priority)*waiting_to_hospital(i);
        
    else
        
        waiting_on_scene(i)=inf;
        waiting_to_hospital(i)=inf;
        waiting_on_scene_penalized(i)=inf;
        waiting_to_hospital_penalized(i)=inf;
        which_ambulance(i)=0;
        
    end
end







