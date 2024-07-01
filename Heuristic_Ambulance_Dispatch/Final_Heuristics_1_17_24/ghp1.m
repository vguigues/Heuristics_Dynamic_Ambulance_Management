

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Vincent Guigues.
%Creation date: 15/7/2020.
%Last updated: 21/7/2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: given a list of calls, allocates ambulances to calls using
%Heuristic GHP1 of the paper "New Heuristics for the Operation of an Ambulance 
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

%ambulances(i).hospital_destination.lat is the latitude of the last hospital
%ambulance i went to or is currently forecast to go.

%ambulances(i).hospital_destination.long is the longitude of the last hospital
%ambulance i went to or is currently forecast to go.

%ambulance(i).arrival_time_at_f_last_trip: time ambulance i was freed from
%the last hospital visited or time it will be freed from an hospital if it is
%currently going to an hospital. If the call is not sent to hospital this
%is the time the ambulance becomes available after serving the call on the scene

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
%5: stays at hospital waiting to be freed, 6: going from hospital to base.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bases(i).lat is the latitude of base i
%bases(i).long is the longitude of base i

%nb_base: number of bases

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

function [waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_noforward_priorities_closest_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases,penalty_matrix,samu_type,inittime)

calls_end=zeros(1,nb_calls);
waiting_on_scene=zeros(1,nb_calls);
waiting_to_hospital=zeros(1,nb_calls);
waiting_on_scene_penalized=zeros(1,nb_calls);
waiting_to_hospital_penalized=zeros(1,nb_calls);
which_ambulance=zeros(1,nb_calls);

queue_size=0;

%index of the last call treated
index_call=1;

%queue contains the queue of calls. We only need to store in
%queue the indexes of the calls (the data of the calls is stored
%in calls).
%queue(i) is the index of i-th call in the queue
%calls are ordered by priority with queue(1) the most prioritary call
queue=[];

current_time=inittime;

%event_call=1 if the event treated is call and 0 otherwise (in which case
%the event is "an ambulance is freed")

future_arrival_time_hospital=[];
for j=1:nb_ambulances
    if (ambulances(j).arrival_time_at_f_last_trip>current_time)
        future_arrival_time_hospital=[future_arrival_time_hospital;ambulances(j).arrival_time_at_f_last_trip];
    end
end

if (length(future_arrival_time_hospital)>0)
    minh=min(future_arrival_time_hospital);
    if (calls(index_call).time<=minh)
        event_call=1;
        current_time=calls(index_call).time;
    else
        event_call=0;
        current_time=minh;
    end
else
    event_call=1;
    current_time=calls(index_call).time;
end

%Current time, updated at time instants corresponding to our 2 types of
%events: 1) calls and 2) ambulance is freed.

%Number of calls already attended
calls_attended=0;

while (calls_attended<nb_calls)
    if (event_call==1)
        %we add the call index to the queue
        queue=[queue;index_call];
        queue_size=queue_size+1;
    end
    
    if (length(queue_size)>=2)
        %We need to update the order of the calls in the queue taking
        %into account the priorities and penalizations.
        %The calls are ordered in decreasing order of the penalized waiting time:
       
        Times=zeros(queue_size,1);
        for k=1:queue_size
            priority_factor=penalization(calls(queue(k)).priority);
            Times(k)=priority_factor*(current_time-calls(queue(k)).time);
        end
        [val,ind]=sort(Times,'descend');
        queue=queue(ind);
    end
    
    queue_aux=[];
    laux=0;
    for k=1:queue_size
        %We treat the calls by order of priority
        %We allocate an ambulance to a call if this is the quickest to
        %reach the call and it is free
        mintimek=inf;
        mintimekp=inf;
        index_ambs=[];
        for j=1:nb_ambulances
                [bool]=can_answer(ambulances(j).type,calls(queue(k)).priority,samu_type);
                if (bool)
                    if (ambulances(j).arrival_time_at_b_last_trip<=current_time)
                        time_this_ambulance=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                    elseif (ambulances(j).arrival_time_at_f_last_trip<=current_time)
                        %Where is the ambulance?
                        %It started at ambulance(j).arrival_time_at_f_last_trip its
                        %travel from [ambulances(j).hospital_destination.lat,ambulances(j).hospital_destination.long]
                        %to its base [ambulances(j).base.lat,ambulances(j).base.long]
                        %Where is it at current_time?
                        [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).hospital_destination.lat,ambulances(j).hospital_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,current_time,ambulances(j).speed);
                        time_this_ambulance=travel_time(curr_lat,curr_long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                    else
                        %The ambulance is not free if and only if
                        %ambulance(j).arrival_time_at_f_last_trip>calls(queue(k)).time
                        %In this case to go to [calls(queue(k)).lat,calls(queue(k)).long] the ambulance will
                        %start from the hospital it is currently headed to
                        %Time to arrive at hospital for j plus time to go from the hospital to the call
                        time_to_hospital=ambulances(j).arrival_time_at_f_last_trip-current_time;
                        time_from_hosp_to_call=travel_time(ambulances(j).hospital_destination.lat,ambulances(j).hospital_destination.long,calls(queue(k)).lat,calls(queue(k)).long,ambulances(j).speed);
                        time_this_ambulance=time_to_hospital+time_from_hosp_to_call;
                    end
                    time_this_ambulancep=penalized_response_time(penalization,time_this_ambulance,penalty_matrix,samu_type,ambulances(j).type,calls(queue(k)).priority,quality_coeff);
                    
                    if (time_this_ambulancep==mintimekp)
                        index_ambs=[index_ambs;j];
                        mintimek=time_this_ambulance;
                    elseif (time_this_ambulancep<mintimekp)
                        index_ambs=[j];
                        mintimekp=time_this_ambulancep;
                        mintimek=time_this_ambulance;
                    end
                end
        end
        
        boolfound=0;
        index_amb=0;
        for i=1:length(index_ambs)
            if (ambulances(index_ambs(i)).arrival_time_at_f_last_trip<=current_time)
                if (boolfound==0)
                    boolfound=1;
                    index_amb=index_ambs(i);
                elseif (ambulances(index_ambs(i)).type>ambulances(index_amb).type)
                    index_amb=index_ambs(i);         
                end
            end
        end        
        if (index_amb>0)
            [ambulances,waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,ambulances_times,ambulances_trips,trip_type,calls_end]=allocate_ambulance(index_amb,queue(k),ambulances,calls,trip_type,ambulances_times,ambulances_trips,mintimek,hospitals,penalization,cleaning_bases,waiting_on_scene,waiting_to_hospital,calls_end,waiting_on_scene_penalized,waiting_to_hospital_penalized,traveltype,basechoice,Rearth);
            which_ambulance(queue(k))=index_amb;
        else
            queue_aux=[queue_aux;queue(k)];
            laux=laux+1;
        end
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


for i=1:nb_calls
    waiting_on_scene_penalized(i)=penalization(calls(i).priority)*waiting_on_scene(i);
    waiting_to_hospital_penalized(i)=penalization(calls(i).priority)*waiting_to_hospital(i);
end




