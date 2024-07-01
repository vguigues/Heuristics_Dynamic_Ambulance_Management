
%ridetype='geodesic', 'street'
%basechoice='closest', 'home', 'optimized'

function [waiting_to_hospital,waiting_on_scene_penalized,waiting_to_hospital_penalized,ambulances_times,ambulances_trips,trip_type,calls_end,ambulances]=allocate_ambulance_to_call_real_traj(ambulances,index_amb,calls,i,ambulances_times,ambulances_trips,trip_type,graph,mintime,hospitals,nbTotalNodes,bases,penalization,cleaning_bases,ridetype,basechoice)

%If index_amb=0 this  means that the call cannot be attended because no
%ambulance can attend calls of the corresponding priority.
%Such situation should not happen in practice.
if (index_amb>0)
    if (ambulances(index_amb).arrival_time_at_b_last_trip<calls(i).time)
        ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},calls(i).time];
        %ambulances_trips{1,index_amb}
        %ambulances(index_amb).base
        %[graph(ambulances(index_amb).base).lat;graph(ambulances(index_amb).base).long]
        ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(ambulances(index_amb).base).lat;graph(ambulances(index_amb).base).long]];
        trip_type{1,index_amb}=[trip_type{1,index_amb},1];
    elseif (ambulances(index_amb).arrival_time_at_f_last_trip<=calls(i).time)
        if (ambulances(index_amb).arrival_time_at_f_last_trip==calls(i).time)
            len=length(ambulances_times{1,index_amb});
            ambulances_trips{1,index_amb}=ambulances_trips{1,index_amb}(:,1:len-1);
            ambulances_times{1,index_amb}=ambulances_times{1,index_amb}(:,1:len-1);
            trip_type{1,index_amb}=trip_type{1,index_amb}(:,1:len-2);
        else
            %path_free_to_base.l path_free_to_base.t path_free_to_base.n
            %path_to_call.l path_to_call.t path_to_call.n
            [type,node1,node2,lat,long,nbTotalNodes]=manageposition(ambulances(index_amb).path_free_to_base.l,ambulances(index_amb).path_free_to_base.t(1),ambulances(index_amb).path_free_to_base.n,calls(i).time,ambulances(index_amb).speed,graph,nbTotalNodes);
            len=length(ambulances_times{1,index_amb});
            ambulances_trips{1,index_amb}(:,len)=[lat;long];
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
    
    if (calls(i).cleaning_needed)
        if (calls(i).hosp_needed)
            [tAux,L]=street_travel_time(calls(i).id,hospitals(calls(i).ih).id,ambulances(index_amb).speed,graph);
            t_from_c_to_h=tAux;
            [tAux,L]=street_travel_time(hospitals(calls(i).ih).id,hospitals(calls(i).ih).cbase,ambulances(index_amb).speed,graph);
            t_from_h_to_cb=tAux;
            [tAux,L]=street_travel_time(hospitals(calls(i).ih).cbase,cleaning_bases(hospitals(calls(i).ih).cbaseindex).base,ambulances(index_amb).speed,graph);
            t_from_cb_to_b=tAux;
            ambulances(index_amb).path_free_to_base.l=L;
            ambulances(index_amb).path_free_to_base.n=length(L);
            tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
            tArrivalcb=tLeaveH+t_from_h_to_cb;
            tLeavecb=tArrivalcb+calls(i).cleaning_time;
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
            ambulances(index_amb).base=cleaning_bases(hospitals(calls(i).ih).cbaseindex).base;
            
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
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(cleaning_bases(hospitals(calls(i).ih).cbaseindex).base).lat,graph(cleaning_bases(hospitals(calls(i).ih).cbaseindex).base).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
        else
            mindist=inf;
            indmin=0;
            for indk=1:length(cleaning_bases)
                [d,listnodes]=dijkstraList(graph,calls(i).id,cleaning_bases(indk).id);
                if (d<mindist)
                    mindist=d;
                    indmin=indk;
                end
            end
            cb=cleaning_bases(indmin).id;
            thisBase=cleaning_bases(indmin).base;
            [t_from_c_to_cb,listnodes]=street_travel_time(calls(i).id,cb,ambulances(index_amb).speed,graph);
            [t_from_cb_to_b,listnodes]=street_travel_time(cb,thisBase,ambulances(index_amb).speed,graph);
            timeArrivalcb=timeLeaveScene+t_from_c_to_cb;
            timeLeavecb=timeArrivalcb+calls(i).cleaning_time;
            timeBase=timeLeavecb+t_from_cb_to_b;
            waiting_to_hospital=calls(i).time_on_scene;
            
            ambulances(index_amb).path_free_to_base.l=listnodes;
            ambulances(index_amb).path_free_to_base.n=length(listnodes);
            ambulances(index_amb).path_free_to_base.t(1)=timeLeavecb;
            for indL=2:length(listnodes)
                whichInd=find(graph(listnodes(indL-1)).neighbors==listnodes(indL));
                ambulances(index_amb).path_free_to_base.t(indL)=ambulances(index_amb).path_free_to_base.t(indL-1)+(graph(listnodes(indL-1)).distances(whichInd)/ambulances(index_amb).speed);
            end
            
            %ambulances(index_amb).free_destination.lat=graph(cb).lat;
            %ambulances(index_amb).free_destination.long=graph(cb).long;
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
            
            calls_end=timeLeaveScene;
        end
    else
        if (calls(i).hosp_needed)
                                                    
            [tAux,L]=street_travel_time(calls(i).id,hospitals(calls(i).ih).id,ambulances(index_amb).speed,graph);
            t_from_c_to_h=tAux;
            [tAux,L]=street_travel_time(hospitals(calls(i).ih).id,hospitals(calls(i).ih).base,ambulances(index_amb).speed,graph);
            t_from_h_to_b=tAux;
                    
            tLeaveH=timeLeaveScene+calls(i).time_at_hospital+t_from_c_to_h;
            tBase=tLeaveH+t_from_h_to_b;
            ambulances(index_amb).path_free_to_base.l=L;
            ambulances(index_amb).path_free_to_base.n=length(L);
            ambulances(index_amb).path_free_to_base.t(1)=tLeaveH;
            
            
            %for indL=2:length(L)
            %    whichInd=find(graph(listnodes(indL-1)).neighbors==L(indL));
            %    ambulances(index_amb).path_free_to_base.t(indL)=ambulances(index_amb).path_free_to_base.t(indL-1)+(graph(L(indL-1)).distances(whichInd)/ambulances(index_amb).speed);
            %end
            for index=2:ambulances(index_amb).path_free_to_base.n
                indA=find(graph(L(index-1)).neighbors==L(index));
                ambulances(index_amb).path_free_to_base.t(index)=ambulances(index_amb).path_free_to_base.t(index-1)+(graph(L(index-1)).distances(indA))/ambulances(index_amb).speed;
            end
            waiting_to_hospital=calls(i).time_on_scene+t_from_c_to_h;
            ambulances(index_amb).free_destination=hospitals(calls(i).ih).id;
            ambulances(index_amb).arrival_time_at_f_last_trip=tLeaveH;
            ambulances(index_amb).arrival_time_at_b_last_trip=tBase;
            ambulances(index_amb).base=hospitals(calls(i).ih).base;
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH-calls(i).time_at_hospital];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},4];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tLeaveH];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[hospitals(calls(i).ih).lat;hospitals(calls(i).ih).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},5];
            
            ambulances_times{1,index_amb}=[ambulances_times{1,index_amb},tBase];
            ambulances_trips{1,index_amb}=[ambulances_trips{1,index_amb},[graph(hospitals(calls(i).ih).base).lat;graph(hospitals(calls(i).ih).base).long]];
            trip_type{1,index_amb}=[trip_type{1,index_amb},8];
            calls_end=timeLeaveScene+t_from_c_to_h;
        else
            mindist=inf;
            indmin=0;
            for indk=1:nb_base
                [d]=street_distance(calls(i).id,bases(indk).id,graph);
                if (d<mindist)
                    mindist=d;
                    indmin=indk;
                end
            end

            [t_from_c_to_b,L]=street_travel_time(calls(i).id,bases(indmin).id,ambulances(index_amb).speed,graph);
            
            timeBase=timeLeaveScene+t_from_c_to_b;
            waiting_to_hospital=calls(i).time_on_scene;
            %ambulances(index_amb).free_destination.lat=calls(i).lat;
            %ambulances(index_amb).free_destination.long=calls(i).long;
            ambulances(index_amb).free_destination=calls(i).id;
            ambulances(index_amb).arrival_time_at_f_last_trip=timeLeaveScene;
            ambulances(index_amb).arrival_time_at_b_last_trip=timeBase;
            %ambulances(index_amb).base.lat=bases(indmin).lat;
            ambulances(index_amb).base=bases(indmin).id;
            
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
            
            calls_end=timeLeaveScene;
        end
    end
    waiting_on_scene_penalized=penalization(calls(i).priority)*mintime;
    waiting_to_hospital_penalized=penalization(calls(i).priority)*waiting_to_hospital;
else    
    waiting_to_hospital=inf;
    waiting_on_scene_penalized=inf;
    waiting_to_hospital_penalized=inf;
end
