
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Vincent Guigues.
%Creation date: 1/10/2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: given a list of calls, allocates ambulances to calls using
%Section 7 of the paper "New Heuristics for the Operation of an Ambulance 
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
%integer, the more the call is prioritary.

%calls(i).time_on_scene: time the ambulance stays on the scene of the call

%calls(i).hosp_needed=1 if the patient of the call needs to go to hospital
%and 0 otherwise.

%calls(i).time_at_hospital: time the ambulance waits at hospital before
%being available for another dispatch.

%calls(i).cleaning_needed=1 if the ambulance in charge of call i needs
%cleaning after service and 0 otherwise.

%calls(i).cleaning_time: time required to clean the ambulance attending
%call i.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ambulances is a list of structures with ambulances(i) a structure for the
%data of ambulance i with the following fields:

%ambulances(i).base.lat is the latitude of the base of ambulance i

%ambulances(i).base.long is the longitude of the base of ambulance i

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
%at its base either for its last trip if it is already at its base or
%the future arrival time at the base if it is currently in service.

%ambulances(i).speed is the speed of the ambulance.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb_ambulances: number of ambulances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hospitals: a list of structures giving the locations of hospitals
%hospitals(i).lat is the the latitude of hospital i
%hospitals(i).long is the the longitude of hospital i
%hospitals(i).cbase.lat is the latitude of the closest clenaing base to hospital i
%hospitals(i).cbase.long is the longitude of the closest base to hospital i

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
%5: stays at hospital waiting to be freed, 6: going to a cleaning base.
%7: stays at the cleaning base. 8: goes back to a parking base.

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


%calls_end(i) is the instant the patient of call i arrives at hospital.

function [waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_forward_nomyopic_closest_basent(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases)

calls_end=zeros(1,nb_calls);
waiting_on_scene=zeros(1,nb_calls);
waiting_to_hospital=zeros(1,nb_calls);
waiting_on_scene_penalized=zeros(1,nb_calls);
waiting_to_hospital_penalized=zeros(1,nb_calls);
which_ambulance=zeros(1,nb_calls);
i=1;
isAllocated=zeros(nb_calls,1);
travelTimes=inf*ones(nb_calls,nb_ambulances);
maxIndex=0;
mintime=inf(nb_calls,1);
index_amb=zeros(nb_calls,1);
amb_type=zeros(nb_calls,1);
inds=[];

  1 2  3   4  5  6    minp
1   D  Db            
2   Db    Db    Db   20  
3   Db               30>20 2 para 3 
4   Db    Db    D        
5         Db          30>20

tc<=tf


consideramos uma melhor ambulancia j 
j=2
j=4
considerar as chamadas 
tc<=tf
tal que j e tbem uma das melhores ambulancias

while (i<=nb_calls)
    
    if (isAllocated(i))
        i=i+1;
    else
        %We will allocate an ambulance for call i and possibly for future
        %calls that would require the same "best ambulance" as call i.
        %Minimal response time for call i
        %Type of best ambulance for call i
        %Compute travel times forward until call index maxIndex
        
        if (i==(maxIndex+1))
            disp('maxIndex+1=i')
            i
            %We allocate an ambulance to call with index i
            for j=1:nb_ambulances
                %First check if ambulance j can take calls of priority calls(i).priority)
                if (ambulances(j).type<=calls(i).priority)
                    if (ambulances(j).arrival_time_at_b_last_trip<=calls(i).time)
                        travelTimes(i,j)=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(i).lat,calls(i).long,ambulances(j).speed);
                    elseif (ambulances(j).arrival_time_at_f_last_trip<=calls(i).time)
                        [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,calls(i).time,ambulances(j).speed);
                        travelTimes(i,j)=travel_time(curr_lat,curr_long,calls(i).lat,calls(i).long,ambulances(j).speed);
                    else
                        time_to_hospital=ambulances(j).arrival_time_at_f_last_trip-calls(i).time;
                        time_from_hosp_to_call=travel_time(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,calls(i).lat,calls(i).long,ambulances(j).speed);
                        travelTimes(i,j)=time_to_hospital+time_from_hosp_to_call; 
                    end
                    %Improvements: store all ambulances that give minimal response time and then check
                    %below if one is available
                    if ((travelTimes(i,j)<mintime(i))||((travelTimes(i,j)==mintime(i))&&(amb_type(i)<ambulances(j).type)))
                        index_amb(i)=j;
                        amb_type(i)=ambulances(j).type;
                        disp('mintime 1');
                        index_amb(i)
                        mintime(i)=travelTimes(i,j);
                        mintime(i)
                    end
                end
            end
            maxIndex=i;
        end
        
        while (isAllocated(i)==0)
            i
            if (ambulances(index_amb(i)).arrival_time_at_f_last_trip<=calls(i).time)
                isAllocated(i)=1;
                i
                disp('before')
                ambulances_times{1,index_amb(i)}
                trip_type{1,index_amb(i)}
                pause
                [ambulances,waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=allocate_ambulance(index_amb(i),i,ambulances,calls,trip_type,ambulances_times,ambulances_trips,mintime(i),hospitals,penalization,cleaning_bases,waiting_on_scene,waiting_to_hospital,calls_end,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance);
                disp('after')
                ambulances_times{1,index_amb(i)}
                trip_type{1,index_amb(i)}
                maxIndex
                i
                pause
                
                for k=(i+1):maxIndex
                    k
                    maxIndex
                    %Update travel_time(k,index_amb(i)), mintime(k) and index_amb(k)
                    if (isAllocated(k)==0)
                        if (ambulances(index_amb(i)).type<=calls(k).priority)
                            if (ambulances(index_amb(i)).arrival_time_at_b_last_trip<=calls(k).time)
                                travelTimes(k,index_amb(i))=travel_time(ambulances(index_amb(i)).base.lat,ambulances(index_amb(i)).base.long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                            elseif (ambulances(index_amb(i)).arrival_time_at_f_last_trip<=calls(k).time)
                                [curr_lat,curr_long]=position_between_origin_destination(ambulances(index_amb(i)).free_destination.lat,ambulances(index_amb(i)).free_destination.long,ambulances(index_amb(i)).base.lat,ambulances(index_amb(i)).base.long,ambulances(index_amb(i)).arrival_time_at_f_last_trip,calls(k).time,ambulances(index_amb(i)).speed);
                                travelTimes(k,index_amb(i))=travel_time(curr_lat,curr_long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                            else
                                time_to_hospital=ambulances(index_amb(i)).arrival_time_at_f_last_trip-calls(k).time;
                                time_from_hosp_to_call=travel_time(ambulances(index_amb(i)).free_destination.lat,ambulances(index_amb(i)).free_destination.long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                travelTimes(k,index_amb(i))=time_to_hospital+time_from_hosp_to_call;
                            end
                            index_amb(k)=0;
                            mintime(k)=inf;
                            amb_type(k)=0;
                            for j=1:nb_ambulances
                                if (ambulances(j).type<=calls(k).priority)
                                    if ((travelTimes(k,j)<mintime(k))||((travelTimes(k,j)==mintime(k))&&(amb_type(k)<ambulances(j).type)))
                                        index_amb(k)=j;
                                        amb_type(k)=ambulances(j).type;
                                        disp('mintime 2');
                                        mintime(k)
                                        mintime(k)=travelTimes(k,j);
                                    end
                                end
                            end
                        end
                    end
                end
            else
                k=maxIndex+1;
                cont=1;
                while (cont)
                    if (k>nb_calls)
                        cont=0;
                    else
                        if (calls(k).time>ambulances(index_amb(i)).arrival_time_at_f_last_trip)
                            cont=0;
                        else
                            if (isAllocated(k)==0)
                                for j=1:nb_ambulances
                                    %First check if ambulance j can take calls of priority calls(i).priority)
                                    if (ambulances(j).type<=calls(k).priority)
                                        if (ambulances(j).arrival_time_at_b_last_trip<=calls(k).time)
                                            travelTimes(k,j)=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(k).lat,calls(k).long,ambulances(j).speed);
                                        elseif (ambulances(j).arrival_time_at_f_last_trip<=calls(k).time)
                                            [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,calls(k).time,ambulances(j).speed);
                                            travelTimes(k,j)=travel_time(curr_lat,curr_long,calls(k).lat,calls(k).long,ambulances(j).speed);
                                        else
                                            time_to_hospital=ambulances(j).arrival_time_at_f_last_trip-calls(k).time;
                                            time_from_hosp_to_call=travel_time(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,calls(k).lat,calls(k).long,ambulances(j).speed);
                                            travelTimes(k,j)=time_to_hospital+time_from_hosp_to_call;
                                        end
                                    end
                                end
                                mintime(k)=inf;
                                amb_type(k)=0;
                                index_amb(k)=0;
                                for j=1:nb_ambulances
                                    %First check if ambulance j can take calls of priority calls(i).priority)
                                    if (ambulances(j).type<=calls(k).priority)
                                        travelTimes(k,j)
                                        if ((travelTimes(k,j)<mintime(k))||((travelTimes(k,j)==mintime(k))&&(amb_type(k)<ambulances(j).type)))
                                            index_amb(k)=j;
                                            j
                                            amb_type(k)=ambulances(j).type;
                                            disp('mintime 3');
                                            mintime(k)=travelTimes(k,j);
                                            mintime(k)
                                        end
                                    end
                                end
                            end
                            k=k+1;
                        end
                    end
                end
                maxIndex=k-1;
                maxmin=-1;
                whichk=0;
                t1=penalization(calls(i).priority)*mintime(i);
                k=i+1;
                if (i<nb_calls)
                    if (calls(k).time<=ambulances(index_amb(i)).arrival_time_at_f_last_trip)
                        contd=1;
                    else
                        contd=0;
                    end
                    while contd
                        if (isAllocated(k)==0)
                            t2=penalization(calls(k).priority)*mintime(k);
                            if ((mintime(k)==travelTimes(k,index_amb(i)))&&(t2>maxmin))
                                whichk=k;
                                maxmin=t2;
                            end
                        end
                        if (k+1>nb_calls)
                            contd=0;
                        elseif (calls(k+1).time>ambulances(index_amb(i)).arrival_time_at_f_last_trip)
                            contd=0;
                        else
                            k=k+1;
                        end
                    end
                end
                if (maxmin>t1)
                    isAllocated(whichk)=1;
                    disp('before whichk')
                    ambulances_times{1,index_amb(i)}
                    trip_type{1,index_amb(i)}
                    whichk
                    pause
                    
                    %Allocate ambulance with index index_amb(i) to call with index whichk
                    [ambulances,waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=allocate_ambulance(index_amb(i),whichk,ambulances,calls,trip_type,ambulances_times,ambulances_trips,mintime(whichk),hospitals,penalization,cleaning_bases,waiting_on_scene,waiting_to_hospital,calls_end,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance);
                    disp('after whichk')
                    ambulances_times{1,index_amb(i)}
                    trip_type{1,index_amb(i)}
                    whichk
                    i
                    maxIndex
                    pause
                    
                    %Update travel times with ambulance with index index_amb(i)
                    for k=i:maxIndex
                        k
                        isAllocated(k)
                        %Update travel_time(k,index_amb(i)), mintime(k) and index_amb(k)
                        if (isAllocated(k)==0)
                            if (ambulances(index_amb(i)).type<=calls(k).priority)
                                if (ambulances(index_amb(i)).arrival_time_at_b_last_trip<=calls(k).time)
                                    travelTimes(k,index_amb(i))=travel_time(ambulances(index_amb(i)).base.lat,ambulances(index_amb(i)).base.long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                elseif (ambulances(index_amb(i)).arrival_time_at_f_last_trip<=calls(k).time)
                                    [curr_lat,curr_long]=position_between_origin_destination(ambulances(index_amb(i)).free_destination.lat,ambulances(index_amb(i)).free_destination.long,ambulances(index_amb(i)).base.lat,ambulances(index_amb(i)).base.long,ambulances(index_amb(i)).arrival_time_at_f_last_trip,calls(k).time,ambulances(index_amb(i)).speed);
                                    travelTimes(k,index_amb(i))=travel_time(curr_lat,curr_long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                else
                                    time_to_hospital=ambulances(index_amb(i)).arrival_time_at_f_last_trip-calls(k).time;
                                    time_from_hosp_to_call=travel_time(ambulances(index_amb(i)).free_destination.lat,ambulances(index_amb(i)).free_destination.long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                    travelTimes(k,index_amb(i))=time_to_hospital+time_from_hosp_to_call;
                                end
                            end
                            index_amb(k)=0;
                            mintime(k)=inf;
                            amb_type(k)=0;
                            for j=1:nb_ambulances
                                if (ambulances(j).type<=calls(k).priority)
                                    if ((travelTimes(k,j)<mintime(k))||((travelTimes(k,j)==mintime(k))&&(amb_type(k)<ambulances(j).type)))
                                        index_amb(k)=j;
                                        amb_type(k)=ambulances(j).type;
                                        disp('mintime 4')
                                        index_amb(i)
                                        j
                                        mintime(k)=travelTimes(k,j);
                                        mintime(k)
                                    end
                                end
                            end
                            
                        end
                    end
                else
                    isAllocated(i)=1;
                    disp('before i')
                    ambulances_times{1,index_amb(i)}
                    trip_type{1,index_amb(i)}
                    i
                    maxIndex
                    calls(i).time
                    mintime(i)
                    calls(i).time+mintime(i)
                    pause
                    [ambulances,waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=allocate_ambulance(index_amb(i),i,ambulances,calls,trip_type,ambulances_times,ambulances_trips,mintime(i),hospitals,penalization,cleaning_bases,waiting_on_scene,waiting_to_hospital,calls_end,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance);
                    disp('after i')
                    ambulances_times{1,index_amb(i)}
                    trip_type{1,index_amb(i)}
                    i
                    pause
                    
                    for k=i:maxIndex
                        %Update travel_time(k,index_amb(i)), mintime(k) and index_amb(k)
                        if (isAllocated(k)==0)
                            if (ambulances(index_amb(i)).type<=calls(k).priority)
                                if (ambulances(index_amb(i)).arrival_time_at_b_last_trip<=calls(k).time)
                                    travelTimes(k,index_amb(i))=travel_time(ambulances(index_amb(i)).base.lat,ambulances(index_amb(i)).base.long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                elseif (ambulances(index_amb(i)).arrival_time_at_f_last_trip<=calls(k).time)
                                    [curr_lat,curr_long]=position_between_origin_destination(ambulances(index_amb(i)).free_destination.lat,ambulances(index_amb(i)).free_destination.long,ambulances(index_amb(i)).base.lat,ambulances(index_amb(i)).base.long,ambulances(index_amb(i)).arrival_time_at_f_last_trip,calls(k).time,ambulances(index_amb(i)).speed);
                                    travelTimes(k,index_amb(i))=travel_time(curr_lat,curr_long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                else
                                    time_to_hospital=ambulances(index_amb(i)).arrival_time_at_f_last_trip-calls(k).time;
                                    time_from_hosp_to_call=travel_time(ambulances(index_amb(i)).free_destination.lat,ambulances(index_amb(i)).free_destination.long,calls(k).lat,calls(k).long,ambulances(index_amb(i)).speed);
                                    travelTimes(k,index_amb(i))=time_to_hospital+time_from_hosp_to_call;
                                end
                            end
                            index_amb(k)=0;
                            mintime(k)=inf;
                            amb_type(k)=0;
                            for j=1:nb_ambulances
                                if (ambulances(j).type<=calls(k).priority)
                                    disp('BEF')
                                    (travelTimes(k,j)<mintime(k))||((travelTimes(k,j)==mintime(k))&&(amb_type(k)<ambulances(j).type))
                                    travelTimes(k,j)
                                    k
                                    mintime(k)
                                    disp('AFTER')
                                    pause
                                    if ((travelTimes(k,j)<mintime(k))||((travelTimes(k,j)==mintime(k))&&(amb_type(k)<ambulances(j).type)))
                                        index_amb(k)=j;
                                        amb_type(k)=ambulances(j).type;
                                        mintime(k)=travelTimes(k,j);
                                        disp('min time 5')
                                        mintime(k)
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end
        end
    end
end







