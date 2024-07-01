Instance data:
	nb_scenarios: # of scenarios
	nb_hospitals: # of hospitals
	nb_bases: # of ambulance bases
	nb_types_ambulance: # of types of ambulance
	nb_ambulances: # of ambulances
	nb_priorities: # of call types
	x_min, x_max, y_min, y_max: coordinates of simulation region
	A(c): set of ambulance types that can answer to a call of type c
	Calls -> call.time, call.(lat,long), call.type: time, location and priority of calls.
	Hospitals -> hospital.(lat,long): location of hospitals
	Bases -> base.(lat,long), base.capacity: location and capacity of ambulance bases.
	Ambulances-> amb.type, amb.(lat,long), amb.speed: Type, initial location and speed of 													ambulances.

Solver data:
	Instance data: see above
	t0: current time in the simulation
	Call: Received call, null if the problem is about returning an ambulance.
	Ambulance: Ambulance that arrived at some hospital, null if the problem
							 	is about responding a call.
	c0: type of the received call, -1 if the problem is about returning an ambulance.
	l0: Discrete location of received call, null if the problem is about returning an 
		ambulance.
	Ambulances: List of ambulances with their current locations in t0
	H(c,l): set of hospitals that can respond to a call of type c coming from discrete
			location l
	C(c,l): # of calls of type c from discrete location l in queue.
	A0(a,b): # of ambulances of type a at base b in t0
	A0(a,h): # of ambulances of type a at hospital h in t0
	A0(a,l,b): # of ambulances of type a at discrete location l returning to base b in t0
	A0(c,a,l,h): # of ambulances of type a at discrete location l responding to call of type c
				heading towards hospital h
	A0(c,a,l1,l,h): # of ambulances of type a at discrete location l1 responding to call of 				type c at location l and heading towards hospital h after reaching l
	L(t,a,b): set of intermediate locations for ambulances of type a by leaving all hospitals
			at time t towards base b.

	lambda(t,c,l): # of calls of type c arriving at discrete time t in discrete location l.


