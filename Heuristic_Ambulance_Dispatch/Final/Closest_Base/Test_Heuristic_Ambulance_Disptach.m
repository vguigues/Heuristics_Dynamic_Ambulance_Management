
%addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Data_Processing'
%addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Final\Closest_Base'

function [waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type]=Test_Heuristic_Ambulance_Disptach(nb_calls,xmax,ymax,Tmax,prob_to_hosp,penalization,heuristic)

nb_calls=4;
xmax=10;
ymax=10;
penalization=[4;2;1];
Tmax=10;
prob_to_hosp=1;

nb_hospitals=3;
hospitals(1).lat=xmax/2;
hospitals(1).long=0;
hospitals(2).lat=xmax/2;
hospitals(2).long=ymax/2;
hospitals(3).lat=xmax/2;
hospitals(3).long=ymax;

cleaning_bases(1).lat=xmax/4;
cleaning_bases(1).long=ymax/4;
cleaning_bases(2).lat=3*xmax/4;
cleaning_bases(2).long=ymax/4;
cleaning_bases(3).lat=xmax/4;
cleaning_bases(3).long=3*ymax/4;
cleaning_bases(4).lat=3*xmax/4;
cleaning_bases(4).long=3*ymax/4;

times=Tmax*rand(nb_calls,1);
times=sort(times);
for i=1:nb_calls
    calls(i).time=times(i);
    calls(i).lat=xmax*rand;
    calls(i).long=ymax*rand;
    indexhosp=1;
    mindist=mydistance(calls(i).lat,calls(i).long,hospitals(1).lat,hospitals(1).long);
    for k=2:nb_hospitals
        dist=mydistance(calls(i).lat,calls(i).long,hospitals(k).lat,hospitals(k).long);
        if (dist<mindist)
            mindist=dist;
            indexhosp=k;
        end
    end
    
    if (mod(i,10)==0)
        calls(i).cleaning_needed=1;
        calls(i).cleaning_time=0.5;       
    else
        calls(i).cleaning_needed=0;
        calls(i).cleaning_time=0;
    end
       
    calls(i).ih=indexhosp;
    calls(i).priority=mod(i,3)+1;
    %calls(i).priority=3;    
    calls(i).time_on_scene=0.2;
    u=rand;
    if (u<=prob_to_hosp)
        calls(i).hosp_needed=1;
        calls(i).time_at_hospital=0.2;
    else
        calls(i).hosp_needed=0;
        calls(i).time_at_hospital=0;
    end
end

%nb_calls=2;
%%calls(1).time=2;
%calls(1).lat=1;
%calls(1).long=9;
%calls(2).time=2.5;
%calls(2).lat=6;
%calls(2).long=10;
%calls(1).priority=3;

%calls(2).priority=3;

nb_base=4;
bases(1).lat=0;
bases(1).long=0;

bases(2).lat=xmax;
bases(2).long=0;

bases(3).lat=0;
bases(3).long=ymax;

bases(4).lat=xmax;
bases(4).long=ymax;

for i=1:nb_hospitals
    
    %Find closest base to hospital(i)
    mindist=inf;
    indmin=0;
    for k=1:nb_base
        [d]=mydistance(hospitals(i).lat,hospitals(i).long,bases(k).lat,bases(k).long);
        if (d<mindist)
            mindist=d;
            indmin=k;
        end
    end    
    hospitals(i).base.lat=bases(indmin).lat;
    hospitals(i).base.long=bases(indmin).long;
    
    %Find closest cleaning base to hospital(i)
    indmin=0;
    mindist=inf;
    for k=1:length(cleaning_bases)
        [d]=mydistance(hospitals(i).lat,hospitals(i).long,cleaning_bases(k).lat,cleaning_bases(k).long);
        if (d<mindist)
            mindist=d;
            indmin=k;
        end
    end    
    hospitals(i).cbase=indmin;
end

for i=1:length(cleaning_bases)
        mindist=inf;
        for k=1:length(bases)
            [d]=mydistance(cleaning_bases(i).lat,cleaning_bases(i).long,bases(k).lat,bases(k).long);
            if (d<mindist)
                mindist=d;
                indmin=k;
            end
        end
        cleaning_bases(i).base.lat=bases(indmin).lat;
        cleaning_bases(i).base.long=bases(indmin).long;
end    
    

%axes('pos',[0 0.5 .05 .05])
%imshow('C:\Users\vince\Dropbox\Cours\Proba\Cours\Simu_Airplane_Monty_Hall\amb_icon.jpg')

nb_base=4;
nb_ambulances=4;
ambulances(1).base.lat=0;
ambulances(1).base.long=0;
ambulances(1).speed=1;

ambulances(2).base.lat=xmax;
ambulances(2).base.long=0;
ambulances(2).speed=1;

ambulances(3).base.lat=0;
ambulances(3).base.long=ymax;
ambulances(3).speed=1;

ambulances(4).base.lat=xmax;
ambulances(4).base.long=ymax;
ambulances(4).speed=1;

ambulances(1).type=1;
ambulances(2).type=2;
ambulances(3).type=3;
ambulances(4).type=2;

%ambulance(1) is going to its base
%The ambulance comes from hospital 2 and left the hospital at t=-3;
ambulances(1).free_destination.lat=hospitals(2).lat; 
ambulances(1).free_destination.long=hospitals(2).long; 
ambulances(1).arrival_time_at_f_last_trip=-3;
%When will it arrive at the base?
[tbase]=travel_time(hospitals(2).lat,hospitals(2).long,ambulances(1).base.lat,ambulances(1).base.long,ambulances(1).speed);
ambulances(1).arrival_time_at_b_last_trip=ambulances(1).arrival_time_at_f_last_trip+tbase;

ambulances_times=cell(1,nb_ambulances);
ambulances_trips=cell(1,nb_ambulances);
trip_type=cell(1,nb_ambulances);

ambulances_times{1,1}=[0,ambulances(1).arrival_time_at_b_last_trip];
[currentlat,currentlong]=position_between_origin_destination(ambulances(1).free_destination.lat,ambulances(1).free_destination.long,ambulances(1).base.lat,ambulances(1).base.lat,ambulances(1).arrival_time_at_f_last_trip,0,ambulances(1).speed);
ambulances_trips{1,1}=[currentlat;currentlong];
ambulances_trips{1,1}=[ambulances_trips{1,1},[ambulances(1).base.lat;ambulances(1).base.long]];
trip_type{1,1}=6;

%We place ambulance(2) at its base at t=0

ambulances(2).free_destination.lat=[]; 
ambulances(2).free_destination.long=[]; 
ambulances(2).arrival_time_at_f_last_trip=-1;
ambulances(2).arrival_time_at_b_last_trip=0;

ambulances_times{1,2}=[0];
ambulances_trips{1,2}=[10;0];
trip_type{1,2}=[];

%We place ambulance(3) at (2,9) at t=0 and it is going to hospital 3
[thosp]=travel_time(2,9,hospitals(3).lat,hospitals(3).long,ambulances(3).speed);
[tbase]=travel_time(hospitals(3).lat,hospitals(3).long,ambulances(3).base.lat,ambulances(3).base.long,ambulances(3).speed);

ambulances(3).free_destination.lat=hospitals(3).lat; 
ambulances(3).free_destination.long=hospitals(3).long; 
ambulances(3).arrival_time_at_f_last_trip=thosp+0.2;
ambulances(3).arrival_time_at_b_last_trip=thosp+0.2+tbase;

ambulances_times{1,3}=[0,thosp,thosp+0.2,thosp+0.2+tbase];
ambulances_trips{1,3}=[2;9];
ambulances_trips{1,3}=[ambulances_trips{1,3},[hospitals(3).lat;hospitals(3).long]];
ambulances_trips{1,3}=[ambulances_trips{1,3},[hospitals(3).lat;hospitals(3).long]];
ambulances_trips{1,3}=[ambulances_trips{1,3},[ambulances(3).base.lat;ambulances(3).base.long]];
trip_type{1,3}=[4,5,6];

ambulances(4).free_destination.lat=[]; 
ambulances(4).free_destination.long=[]; 
ambulances(4).arrival_time_at_f_last_trip=-1;
ambulances(4).arrival_time_at_b_last_trip=0;

ambulances_times{1,4}=[0];
ambulances_trips{1,4}=[xmax;ymax];
trip_type{1,4}=[];

ambulances_timesp=ambulances_times;
ambulances_tripsp=ambulances_trips;
trip_typep=trip_type;

[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_forward_nomyopic(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases)


%if (heuristic=='fnp')

[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_minmax_priorities_closest_available_basent(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases);

[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_forward_fixed_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base,cleaning_bases);
Plots_Ambulance_Trajectories(nb_ambulances,ambulances_times,ambulances_trips,trip_type,calls,nb_calls,calls_end,ambulances,hospitals,nb_hospitals,bases,nb_base);
%elseif (heuristic=='nfp') 
[waiting_on_scenep,waiting_to_hospitalp,waiting_on_scene_penalizedp,waiting_to_hospital_penalizedp,which_ambulancep,ambulances_timesp,ambulances_tripsp,trip_typep,calls_endp]=closest_among_all_noforward_priorities_fixed_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_timesp,ambulances_tripsp,trip_typep,bases,nb_base);
Plots_Ambulance_Trajectories(nb_ambulances,ambulances_timesp,ambulances_tripsp,trip_typep,calls,nb_calls,calls_endp,ambulances,hospitals,nb_hospitals,bases,nb_base);
%end
[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_available_fixed_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base);

%Closest base
[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_forward_closest_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base);
[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_among_all_noforward_priorities_closest_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base);
[waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,which_ambulance,ambulances_times,ambulances_trips,trip_type,calls_end]=closest_available_closest_base(calls,ambulances,nb_ambulances,hospitals,nb_hospitals,nb_calls,penalization,ambulances_times,ambulances_trips,trip_type,bases,nb_base);

