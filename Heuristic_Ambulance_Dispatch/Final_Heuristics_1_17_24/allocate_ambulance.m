
%Allocates ambulance with index index_amb to call with index i
%mintime is the waiting time between the call arrival and the
%arrival of the ambulance on the scene of the call
%traveltype='G' geodesic rides, traveltype='S' rides along the streets
%basechoice='H' home base, basechoice='C' closest base, basechoice='B' "best" base (heuristic)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%With geodesic rides: 
%ambulances(index_amb).base.lat ambulances(index_amb).base.long
%indexcb=hospitals(i).cbase is the index of the closest cleaning base to
%hospital with index i: this cleaning base is cleaning_bases(indexcb)  
%located at cleaning_bases(indexcb).lat, cleaning_bases(indexcb).long  
%hospitals(i).lat, hospitals(i).long
%calls(i).ih index of the hospital for call i
%ambulances(index_amb).free_destination.lat
%ambulances(index_amb).free_destination.long
%ambulances(index_amb).arrival_time_at_f_last_trip
%ambulances(index_amb).arrival_time_at_b_last_trip
%ambulances(index_amb).base.lat
%ambulances(index_amb).base.long

%cleaning_bases(i).base.lat, cleaning_bases(i).base.long location of the
%base which is the closest to cleaning base with index i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%With rides in streets
%ambulances(index_amb).base is the node id of the base,
%graph(ambulances(index_amb).base).lat,
%graph(ambulances(index_amb).base).long
%ambulances(index_amb).path_free_to_base.l
%ambulances(index_amb).path_free_to_base.t
%ambulances(index_amb).path_free_to_base.n

%ambulances(index_amb).arrival_time_at_f_last_trip
%ambulances(index_amb).arrival_time_at_b_last_trip

%ambulances(index_amb).free_destination is the id of the node
%of the graph corresponding to free_destination

%ambulances(index_amb).base is the id of the node of the graph
%corresponding to base

%hospitals(i).cbase is the id of the node of the graph of the cleaning base
%which is the closest to hospital with index i

%hospitals(i).cbaseindex is the index of the cleaning base closest to
%hospital with index i

%cleaning_bases(i).base is the id of the node of the graph corresponding to
%the base which is the closest to cleaning base with index i

%cleaning_bases(i).id is the id of the node of the graph corresponding to
%cleaning base with index i

%hospitals(i).id is the id of the node of the graph corresponding to
%hospital with index i

%hospitals(i).base is the index of the node of the graph corresponding to
%the base which is the closest to hospital with index i

%calls(i).id is the id of the node of the graph corresponding to call with
%index i

%bases(i).id is the index of the node of the graph corresponding to base
%with index i 

function [ambulances,waiting_on_scene,waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,ambulances_times,ambulances_trips,trip_type,calls_end]=allocate_ambulance(index_amb,i,ambulances,calls,trip_type,ambulances_times,ambulances_trips,mintime,hospitals,penalization,cleaning_bases,waiting_on_scene,waiting_to_hospital,calls_end,waiting_on_scene_penalized,waiting_to_hospital_penalized,traveltype,basechoice,Rearth)

if (ambulances(index_amb).arrival_time_at_b_last_trip<calls(i).time)
    ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},calls(i).time];
    if (traveltype=='G')
        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[ambulances(index_amb).base.lat;ambulances(index_amb).base.long]];
    else
        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(ambulances(index_amb).base).lat;graph(ambulances(index_amb).base).long]];
    end
    trip_type{1,index_amb}=[trip_type{1,index_amb},1];
elseif (ambulances(index_amb).arrival_time_at_f_last_trip<=calls(i).time)
    if  (ambulances(index_amb).arrival_time_at_f_last_trip==calls(i).time)
        len=length(ambulances_times{1,index_amb});
        ambulances_trips{1,index_amb}=ambulances_trips{1,index_amb}(:,1:len-1);
        ambulances_times{1,index_amb}=ambulances_times{1,index_amb}(:,1:len-1);
        trip_type{1,index_amb}=trip_type{1,index_amb}(:,1:len-2);
    else
        len=length(ambulances_times{1,index_amb});
        if (traveltype=='G')
            [curr_lat,curr_long]=position_between_origin_destination_geodesic(ambulances(index_amb).free_destination.lat,ambulances(index_amb).free_destination.long,ambulances(index_amb).base.lat,ambulances(index_amb).base.long,ambulances(index_amb).arrival_time_at_f_last_trip,calls(i).time,ambulances(index_amb).speed,Rearth,'sphere');
            ambulances_trips{1,index_amb}(:,len)=[curr_lat;curr_long];
            ambulances_times{1,index_amb}(len)=calls(i).time; 
        else
            [type,node1,node2,lat,long,nbTotalNodes]=manageposition(ambulances(index_amb).path_free_to_base.l,ambulances(index_amb).path_free_to_base.t(1),ambulances(index_amb).path_free_to_base.n,calls(i).time,ambulances(index_amb).speed,graph,nbTotalNodes,Rearth);
            ambulances_trips{1,index_amb}(:,len)=[lat;long];
            ambulances_times{1,index_amb}(len)=calls(i).time;
        end
    end
else
    len=length(ambulances_times{1,index_amb});
    ambulances_trips{1,index_amb}=ambulances_trips{1,index_amb}(:,1:len-1);
    ambulances_times{1,index_amb}=ambulances_times{1,index_amb}(:,1:len-1);
    trip_type{1,index_amb}=trip_type{1,index_amb}(:,1:len-2);
end

timeArrivalScene=calls(i).time+mintime;
ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeArrivalScene];
if (traveltype=='G')
    ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[calls(i).lat;calls(i).long]];
else
    ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(calls(i).id).lat;graph(calls(i).id).long]];
end
trip_type{1,index_amb}=[trip_type{1,index_amb},2];

timeLeaveScene=timeArrivalScene+calls(i).time_on_scene;
ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeLeaveScene];
if (traveltype=='G')
    ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[calls(i).lat;calls(i).long]];
else
    ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(calls(i).id).lat;graph(calls(i).id).long]];
end
trip_type{1,index_amb}=[trip_type{1,index_amb},3];

waiting_on_scene(i)=mintime;

if (calls(i).cleaning_needed)
    if (calls(i).hosp_needed)
        if (traveltype=='G')
            t_from_c_to_h=travel_time_geodesic(calls(i).lat,calls(i).long,hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,ambulances(index_amb).speed,'sphere',Rearth);
            indexcb=hospitals(calls(i).ih).cbase;
            cbase=cleaning_bases(indexcb);
            t_from_h_to_cb=travel_time_geodesic(hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,cbase.lat,cbase.long,ambulances(index_amb).speed,'sphere',Rearth);
            %Choice of base
            tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
            tArrivalcb=tLeaveH+t_from_h_to_cb;
            tLeavecb=tArrivalcb+calls(i).cleaning_time;
            if (basechoice=='H')
                %We leave cleaning base (cbase.lat,cbase.long) at tLeavecb
                thisbase=ambulances(index_amb).base;
                t_from_cb_to_b=travel_time_geodesic(cbase.lat,cbase.long,thisbase.lat,thisbase.long,ambulances(index_amb).speed,'sphere',Rearth);
                tBase=tLeavecb+t_from_cb_to_b;
            elseif (basechoice=='C')
                thisbase=cleaning_bases(indexcb).base;
                t_from_cb_to_b=travel_time_geodesic(cbase.lat,cbase.long,thisbase.lat,thisbase.long,ambulances(index_amb).speed,'sphere',Rearth);
                tBase=tLeavecb+t_from_cb_to_b;
            else
                [indexBase]=findBestBasegeodesic(ambulances,cbase.lat,cbase.long,tLeavecb,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,Rearth);
                thisbase=bases(indexBase);
                t_from_cb_to_b=travel_time_geodesic(cbase.lat,cbase.long,thisbase.lat,thisbase.long,ambulances(index_amb).speed,'sphere',Rearth);
                tBase=tLeavecb+t_from_cb_to_b;
            end
            
            calls_end(i)=timeLeaveScene+t_from_c_to_h;
            waiting_to_hospital(i)=calls(i).time_on_scene+t_from_c_to_h;
            ambulances(index_amb).free_destination.lat=cbase.lat;
            ambulances(index_amb).free_destination.long=cbase.long;
            ambulances(index_amb).arrival_time_at_f_last_trip=tLeavecb;
            ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
            ambulances(index_amb).base.lat=thisbase.lat;
            ambulances(index_amb).base.long=thisbase.long;
            
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
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[thisbase.lat;thisbase.long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
        else
            [tAux,L]=street_travel_time(calls(i).id,hospitals(calls(i).ih).id,ambulances(index_amb).speed,graph);
            t_from_c_to_h=tAux;
            [tAux,L]=street_travel_time(hospitals(calls(i).ih).id,hospitals(calls(i).ih).cbase,ambulances(index_amb).speed,graph);
            t_from_h_to_cb=tAux;
            tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
            tArrivalcb=tLeaveH+t_from_h_to_cb;
            tLeavecb=tArrivalcb+calls(i).cleaning_time;
            if (basechoice=='H')
               thisbase=ambulances(index_amb).base;
               [tAux,L]=street_travel_time(hospitals(calls(i).ih).cbase,thisbase,ambulances(index_amb).speed,graph);             
            elseif (basechoice=='C')
               thisbase=cleaning_bases(hospitals(calls(i).ih).cbaseindex).base;
               [tAux,L]=street_travel_time(hospitals(calls(i).ih).cbase,thisbase,ambulances(index_amb).speed,graph);
            else 
               [indexBase]=findBestBasestreets(ambulances,hospitals(calls(i).ih).cbase,tLeavecb,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,graph);
               thisbase=bases(indexBase).id;
               [tAux,L]=street_travel_time(hospitals(calls(i).ih).cbase,thisbase,ambulances(index_amb).speed,graph);
            end
            t_from_cb_to_b=tAux;
            ambulances(index_amb).path_free_to_base.l=L;
            ambulances(index_amb).path_free_to_base.n=length(L);
            ambulances(index_amb).path_free_to_base.t(1)=tLeavecb;
            tBase=tLeavecb+t_from_cb_to_b; 
            calls_end=timeLeaveScene+t_from_c_to_h;
            waiting_to_hospital=calls(i).time_on_scene+t_from_c_to_h;
            
            for index=2:ambulances(index_amb).path_free_to_base.n
                indA=find(graph(L(index-1)).neighbors==L(index));
                ambulances(index_amb).path_free_to_base.t(index)=ambulances(index_amb).path_free_to_base.t(index-1)+(graph(L(index-1)).distances(indA))/ambulances(index_amb).speed;
            end
            
            ambulances(index_amb).free_destination=hospitals(calls(i).ih).cbase;
            ambulances(index_amb).arrival_time_at_f_last_trip=tLeavecb;
            ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
            ambulances(index_amb).base=thisbase;
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(i).time_at_hospital];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},4];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},5];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tArrivalcb];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(hospitals(calls(i).ih).cbase).lat;graph(hospitals(calls(i).ih).cbase).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},6];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeavecb];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(hospitals(calls(i).ih).cbase).lat;graph(hospitals(calls(i).ih).cbase).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},7];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(thisbase).lat,graph(thisbase).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            
        end 
    else
        if (traveltype=='G')
            mindist=inf;
            indmin=0;
            for indk=1:length(cleaning_bases)
                [d]=distance_geodesic(calls(i).lat,calls(i).long,cleaning_bases(indk).lat,cleaning_bases(indk).long,'sphere',Rearth)
                if (d<mindist)
                    mindist=d;
                    indmin=indk;
                end
            end
            cb=cleaning_bases(indmin);
            t_from_c_to_cb=travel_time_geodesic(calls(i).lat,calls(i).long,cb.lat,cb.long,ambulances(index_amb).speed,'sphere',Rearth);
            timeArrivalcb=timeLeaveScene+t_from_c_to_cb;
            timeLeavecb=timeArrivalcb+calls(i).cleaning_time;
            if (basechoice=='H')
                    thisbase=ambulances(index_amb).base;
            elseif (basechoice=='C')
                    thisBase=cleaning_bases(indmin).base;    
            else
                    [indexBase]=findBestBasegeodesic(ambulances,cb.lat,cb.long,timeLeavecb,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,Rearth);
                    thisbase=bases(indexBase);
            end
            t_from_cb_to_b=travel_time_geodesic(cb.lat,cb.long,thisBase.lat,thisBase.long,ambulances(index_amb).speed,'sphere',Rearth);
            timeBase=timeLeavecb+t_from_cb_to_b;
            waiting_to_hospital(i)=calls(i).time_on_scene;
            calls_end(i)=timeLeaveScene;
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
        else
            mindist=inf;
            indmin=0;
            for indk=1:length(cleaning_bases)
                [d]=street_distance(calls(i).id,cleaning_bases(indk).id,graph);
                if (d<mindist)
                    mindist=d;
                    indmin=indk;
                end
            end
            cb=cleaning_bases(indmin).id;
            [t_from_c_to_cb,listnodes]=street_travel_time(calls(i).id,cb,ambulances(index_amb).speed,graph);
            timeArrivalcb=timeLeaveScene+t_from_c_to_cb;
            timeLeavecb=timeArrivalcb+calls(i).cleaning_time;
            if (basechoice=='H')
                thisBase=ambulances(index_amb).base
            elseif (basechoice=='C')
                thisBase=cleaning_bases(indmin).base;
            else
               [indexBase]=findBestBasestreets(ambulances,cb,tLeavecb,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,graph);
               thisBase=bases(indexBase).id;
            end
            [t_from_cb_to_b,listnodes]=street_travel_time(cb,thisBase,ambulances(index_amb).speed,graph);

            timeBase=timeLeavecb+t_from_cb_to_b;
            waiting_to_hospital=calls(i).time_on_scene;
            
            ambulances(index_amb).path_free_to_base.l=listnodes;
            ambulances(index_amb).path_free_to_base.n=length(listnodes);
            ambulances(index_amb).path_free_to_base.t(1)=timeLeavecb;
            for indL=2:length(listnodes)
                whichInd=find(graph(listnodes(indL-1)).neighbors==listnodes(indL));
                ambulances(index_amb).path_free_to_base.t(indL)=ambulances(index_amb).path_free_to_base.t(indL-1)+(graph(listnodes(indL-1)).distances(whichInd)/ambulances(index_amb).speed);
            end
            
            ambulances(index_amb).free_destination=cb;
            ambulances(index_amb).arrival_time_at_f_last_trip=timeLeavecb;
            ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
            ambulances(index_amb).base=thisBase;
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeArrivalcb];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(cb).lat;graph(cb).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},6];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeLeavecb];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(cb).lat;graph(cb).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},7];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(thisBase).lat;graph(thisBase).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            
            calls_end(i)=timeLeaveScene;            
        end
    end
else
    if (calls(i).hosp_needed)
        if (traveltype=='G')
            t_from_c_to_h=travel_time_geodesic(calls(i).lat,calls(i).long,hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,ambulances(index_amb).speed,'sphere',Rearth);
            tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
            if (basechoice=='H')
                thisbase=ambulances(index_amb).base;
                t_from_h_to_b=travel_time_geodesic(hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,thisbase.lat,thisbase.long,ambulances(index_amb).speed,'sphere',Rearth);    
            elseif (basechoice=='C')
                thisbase=hospitals(calls(i).ih).base;
                t_from_h_to_b=travel_time_geodesic(hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,thisbase.lat,thisbase.long,ambulances(index_amb).speed,'sphere',Rearth);
            else
                [indexBase]=findBestBasegeodesic(ambulances,hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,tLeaveH,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,Rearth);
                thisbase=bases(indexBase);
                t_from_h_to_b=travel_time_geodesic(hospitals(calls(i).ih).lat,hospitals(calls(i).ih).long,thisbase.lat,thisbase.long,ambulances(index_amb).speed,'sphere',Rearth);
            end            
            tBase=tLeaveH+t_from_h_to_b;
            
            waiting_to_hospital(i)=calls(i).time_on_scene+t_from_c_to_h;
            
            ambulances(index_amb).free_destination.lat=hospitals(calls(i).ih).lat;
            ambulances(index_amb).free_destination.long=hospitals(calls(i).ih).long;
            ambulances(index_amb).arrival_time_at_f_last_trip=tLeaveH;
            ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
            ambulances(index_amb).base.lat=thisbase.lat;
            ambulances(index_amb).base.long=thisbase.long;
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(i).time_at_hospital];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},4];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},5];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[thisbase.lat;thisbase.long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            calls_end(i)=timeLeaveScene+t_from_c_to_h;
        else
            [tAux,L]=street_travel_time(calls(i).id,hospitals(calls(i).ih).id,ambulances(index_amb).speed,graph);
            t_from_c_to_h=tAux;
            tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
            if (basechoice=='H')
               thisbase=ambulances(index_amb).base;
            elseif (basechoice=='C')
               thisbase=hospitals(calls(i).ih).base;
            else 
               [indexBase]=findBestBasestreets(ambulances,hospitals(calls(i).ih).id,tLeaveH,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,graph);
               thisbase=bases(indexBase).id;
            end
            
            [tAux,L]=street_travel_time(hospitals(calls(i).ih).id,thisbase,ambulances(index_amb).speed,graph);
            t_from_h_to_b=tAux;
                    
            tBase=tLeaveH+t_from_h_to_b;
            
            ambulances(index_amb).path_free_to_base.l=L;
            ambulances(index_amb).path_free_to_base.n=length(L);
            ambulances(index_amb).path_free_to_base.t(1)=tLeaveH;
                        
            for index=2:ambulances(index_amb).path_free_to_base.n
                indA=find(graph(L(index-1)).neighbors==L(index));
                ambulances(index_amb).path_free_to_base.t(index)=ambulances(index_amb).path_free_to_base.t(index-1)+(graph(L(index-1)).distances(indA))/ambulances(index_amb).speed;
            end
            
            waiting_to_hospital=calls(i).time_on_scene+t_from_c_to_h;
            
            ambulances(index_amb).free_destination=hospitals(calls(i).ih).id;
            ambulances(index_amb).arrival_time_at_f_last_trip=tLeaveH;
            ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
            ambulances(index_amb).base=thisbase;
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(i).time_at_hospital];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},4];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},5];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(thisbase).lat;graph(thisbase).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            calls_end(i)=timeLeaveScene+t_from_c_to_h;
        end
    else
        %Checked until here no clean no hosp
        if (traveltype=='G')
            if (basechoice=='H')
                thisbase=ambulances(index_amb).base;
                baselat=thisbase.lat;
                baselong=thisbase.long;
                t_from_c_to_b=travel_time_geodesic(calls(i).lat,calls(i).long,baselat,baselong,ambulances(index_amb).speed,'sphere',Rearth);    
            elseif (basechoice=='C')
                mindist=inf;
                indmin=0;
                for indk=1:nb_base
                    [d]=mydistance(calls(i).lat,calls(i).long,bases(indk).lat,bases(indk).long);
                    if (d<mindist)
                        mindist=d;
                        indmin=indk;
                    end
                end
                baselat=bases(indmin).lat;
                baselong=bases(indmin).long;
                t_from_c_to_b=travel_time_geodesic(calls(i).lat,calls(i).long,baselat,baselong,ambulances(index_amb).speed,'sphere',Rearth);    
            else
                [indexBase]=findBestBasegeodesic(ambulances,calls(i).lat,calls(i).long,timeLeaveScene,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,Rearth);
                thisbase=bases(indexBase);
                baselat=thisbase.lat;
                baselong=thisbase.long;
                t_from_c_to_b=travel_time_geodesic(calls(i).lat,calls(i).long,baselat,baselong,ambulances(index_amb).speed,'sphere',Rearth);
            end
            timeBase=timeLeaveScene+t_from_c_to_b;
            waiting_to_hospital(i)=calls(i).time_on_scene;
            ambulances(index_amb).free_destination.lat=calls(i).lat;
            ambulances(index_amb).free_destination.long=calls(i).long;
            ambulances(index_amb).arrival_time_at_f_last_trip=timeLeaveScene;
            ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
            ambulances(index_amb).base.lat=baselat;
            ambulances(index_amb).base.long=baselong;
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[baselat;baselong]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            calls_end(i)=timeLeaveScene;
        else
            if (basechoice=='H')
                thisbase=ambulances(index_amb).base;
                [tAux,L]=street_travel_time(calls(i).id,thisbase,ambulances(index_amb).speed,graph);             
            elseif (basechoice=='C')
                mindist=inf;
                indmin=0;
                for indk=1:nb_base
                     mindist=inf;
                     indmin=0;
                     [d]=street_distance(calls(i).id,bases(indk).id,graph);
                     if (d<mindist)
                         mindist=d;
                         indmin=indk;
                     end
                end     
                thisbase=bases(indmin).id;
                [tAux,L]=street_travel_time(calls(i).id,thisbase,ambulances(index_amb).speed,graph);
            else
               [indexBase]=findBestBasestreets(ambulances,calls(i).id,timeLeaveScene,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,graph);
               thisBase=bases(indexBase).id;
               [tAux,L]=street_travel_time(calls(i).id,thisbase,ambulances(index_amb).speed,graph);
            end
            t_from_c_to_b=tAux; 
            timeBase=timeLeaveScene+t_from_c_to_b;
            waiting_to_hospital=calls(i).time_on_scene;
            %ambulances(index_amb).free_destination.lat=calls(i).lat;
            %ambulances(index_amb).free_destination.long=calls(i).long;
            ambulances(index_amb).free_destination=calls(i).id;
            ambulances(index_amb).arrival_time_at_f_last_trip=timeLeaveScene;
            ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
            ambulances(index_amb).base=thisbase;
            
            ambulances(index_amb).path_free_to_base.l=L;
            ambulances(index_amb).path_free_to_base.n=length(L);
            ambulances(index_amb).path_free_to_base.t(1)=tLeaveH;
            for indL=2:length(L)
                whichInd=find(graph(L(indL-1)).neighbors==L(indL));
                ambulances(index_amb).path_free_to_base.t(indL)=ambulances(index_amb).path_free_to_base.t(indL-1)+(graph(L(indL-1)).distances(whichInd)/ambulances(index_amb).speed);
            end
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},timeBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(bases(indmin).id).lat;graph(bases(indmin).id).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            
            calls_end(i)=timeLeaveScene;
        end
    end
end

waiting_on_scene_penalized(i)=penalization(calls(i).priority)*waiting_on_scene(i);
waiting_to_hospital_penalized(i)=penalization(calls(i).priority)*waiting_to_hospital(i);
     