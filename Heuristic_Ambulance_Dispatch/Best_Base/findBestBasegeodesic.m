

%The ambulance is at (amblat, amblong) at time t0
%The intensity of the Poisson distribution in Voronoi cell with index b 
%and time slot t priority p is lambdaBase(b,t,p)
%The time slot limits are  times(1)  timeArrival times(2) timeArrival+timeAhead times(3)  times(T)
%ambulances(j).base
%ambulances(j).arrival_time_at_b_last_trip

function [indexBase]=findBestBasegeodesic(ambulances,amblat,amblong,t0,lambdaBase,T,nb_base,bases,indexAmb,timeAhead,times,nb_ambulances,poissquant,Rearth)

times=[times,times+times(T)];

for ibase=1:nb_base
    [t]=travel_time_geodesic(amblat,amblong,bases(ibase).lat,bases(ibase).long,ambulances(indexAmb).speed,'sphere',Rearth);
    timeArrival(ibase)=t0+t;
end

%We will test priorities from calls(i).priority to 3
thispriority=ambulances(indexAmb).type;
contd=1;
deficitGlobalMax=-inf;
mintimeG=inf;

while ((thispriority<=3)&&(contd))
    %Computation of maximal ambulance deficit for calls of priority
    %thispriority
    %Candidate ambulances have indexes j such that ambulances(j).type<=thispriority
    deficitMax=-inf;
    mintime=inf;
    for ibase=1:nb_base
        nbAmb=0;
        for index2=1:nb_ambulances
            %Deficit of ambulances at base with index ibase for calls of priority
            %thispriority on time window
            %[timeArrival(ibase),timeArrival(ibase)+timeAhead]
            %How many ambulances will arrive before timeArrival(ibase) at base ibase?
            if (index2~=indexAmb)
                if ((ambulances(index2).arrival_time_at_b_last_trip<=timeArrival(ibase))&&(ambulances(index2).base.lat==bases(ibase).lat)&&(ambulances(index2).base.long==bases(ibase).long)&&(ambulances(index2).type<=thispriority))
                    nbAmb=nbAmb+1;
                end
            end
        end
        index=1;
        while (times(index)<=timeArrival(ibase))
            index=index+1;
        end
        if (times(index)>=(timeArrival(ibase)+timeAhead))
            lambda=timeAhead*lambdaBase(ibase,index-1,thispriority);
        else
            lambda=(times(index)-timeArrival(ibase))*lambdaBase(ibase,index-1,thispriority);
            while (times(index+1)<=timeArrival(ibase)+timeAhead)
                lambda=lambda+(times(index+1)-times(index))*lambdaBase(ibase,index,thispriority);
                index=index+1;
            end
            lambda=(timeArrival(ibase)+timeAhead-times(index))*lambdaBase(ibase,index,thispriority);
        end
        callsquant=poissinv(poissquant,lambda);
        if ((callsquant-nbAmb>deficitMax)||((((callsquant-nbAmb)==deficitMax))&&(timeArrival(ibase)<mintime)))
            deficitMax=callsquant-nbAmb;
            mintime=timeArrival(ibase);
            indexBase=ibase;
        end
        if ((callsquant-nbAmb>deficitGlobalMax)||(((callsquant-nbAmb)==deficitGlobalMax))&&(timeArrival(ibase)<mintimeG))
            mintimeG=timeArrival(ibase);
            deficitGlobalMax=callsquant-nbAmb;
            indexBaseG=ibase;
        end
    end
    if (deficitMax>0)
        contd=0;
    else
        thispriority=thispriority+1;
    end
end

if (deficitGlobalMax<=0)
    indexBase=indexBaseG;
end

