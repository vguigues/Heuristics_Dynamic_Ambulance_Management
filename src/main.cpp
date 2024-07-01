#include "../include/travel.h"
#include "../include/instance.h"
#include "../include/osrm_helper.h"
#include "../include/call_model.h"
#include "../include/policy_tester.h"
#include "../include/solver.h"
#include "../include/data.h"
#include "../include/main.h"

//ILOSTLBEGIN
namespace po = boost::program_options;
std::string config_file;
Param g_params;

class PolicyTester;


int main(int argc, char* argv[]){
	srand(time(NULL));
	
	g_params = load_params(argc,argv);

	std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);


	// GRBEnv env;
	// GeneratorRegressor gen(env);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_old = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("Run time = {}\n", run_time_old);

	// for(int nb = 8; nb <= 30; nb += 2){
	// 	g_params.n_ambulances = nb;
	// 	g_params.n_bases = nb;
	// 	Instance ins;
	// 	Travel travel;
	// 	try{
	// 		PolicyTester pt(ins);
	// 		pt.run();
	// 	}catch(GRBException& ex){
	// 		fmt::print("ERROR {}: {}\n", ex.getErrorCode(), ex.getMessage());
	// 	}

	// 	cin.get();
	// }

	Instance ins;
	Travel travel;
	try{
		PolicyTester pt(ins);
		pt.run();
	}catch(GRBException& ex){
		fmt::print("ERROR {}: {}\n", ex.getErrorCode(), ex.getMessage());
	}


	// Instance ins("instances/no_myopic7.4.txt");
	// ins.travel.euclidian = true;
	// PolicyTester pt(ins);
	// pt.run();


	// Tests Talarico Model
	// GRBEnv env;
	// Instance ins("instances/cont_10_4_4_4.txt");
	// Travel travel(true);
	// // for(auto& amb: ins.ambulances){
	// // 	amb.speed = 1;
	// // }
	// std::shared_ptr<Solver> solver = static_pointer_cast<Solver>(make_shared<GenMinMaxSolver>(
	// 	env,ins.calls[0],ins.ambulances,ins, travel));
	// solver->time = ins.calls[0][8].time+1;
	// fmt::print("Solver time {}\n", solver->time);
	// for(int i = 0; i < 8; ++i){
	// 	std::cout << solver->calls[i] << "\n";
	// 	solver->queue.push_back(i);
	// }
	
	// CallModel cm(env, *solver);
	// cm.solve();
	// cm.print_results();
}


void version(){
	std::cout << "Simulator of ambulance dispatch to emergency calls.\n" <<
		"Version number 0.1 - Developed by: Victor Hugo (UFRJ-PESC-LABOTIM).\n";
}

Param load_params(int argc, char** argv){
	// Declare a group of options that will be 
	// allowed only on command line
	po::options_description generic("Generic Options");
	generic.add_options()
    	("version,v", "Prints version message.")
    	("help,h", "Prints help message.")
    	("file,f", po::value<std::string>(&config_file)->default_value("test.cfg"),
                  "configuration file path.")
    ;
    // Declare a group of options that will be 
	// allowed both on command line and in
	// config file
	po::options_description config("Configuration");
	config.add_options()
	    ("instance,i", po::value<std::string>(), 
	    	"instance file path.")
	    ("solver,s", po::value<std::string>()->default_value("queue"),
	    	"solver used in each call. queue|forward|priorities|minmax|non_miopyc")
	    ("generator_folder,G", po::value<std::string>()->default_value("calibration/Rect10x10"),
	    	"scenario generator file path.")
		("instance_type", po::value<std::string>()->default_value("one_priority"),
			"Type of instance generated.")
		("amb_setup", po::value<std::string>()->default_value("rj"),
			"Ambulance setup. us = 2 ambulance types, 4 priorities. rj = 3 ambulance types, 3 priorities.")
		("closest_base", po::bool_switch()->default_value(true),
			"Whether ambulances should return to closest base or to the best base given deficits. default true")
		("best_base", po::bool_switch()->default_value(false),
			"Whether ambulances should return to best base given statistical data. Overrides closest_base parameter. default false")
	    ("debug,d", po::value<bool>()->default_value(false),
	    	"Runs in debug mode. default = false")
	    ("h_use_fixed_bases,B", po::value<bool>()->default_value(false),
	    	"Ambulance must return to its fixed base. default = false")
	    ("h_forward, F", po::value<bool>()->default_value(false),
	    	"Calls can be dispatched to busy ambulances. default = false")
	    ("h_discard, D", po::value<bool>()->default_value(false),
	    	"Heuristic don't consider a call queue. default = false")
	    ("h_order_priorities, P", po::value<bool>()->default_value(false),
	    	"In Heuristic, order calls in queue by priority. default = false")
	    ("h_order_time, T", po::value<bool>()->default_value(false), 
	    	"In Heuristic, order calls by time. default = false")
	    ("n_nearest_hospitals", po::value<int>()->default_value(1),
	    	"Number of nearest hospitals to a location l that can receive the patients of an emergency in l. default = 1")
	    ("n_nearest_ambs", po::value<int>()->default_value(5),
	    	"In the Model Heuristic, number of nearest ambulances that will be evaluated via optimization model. default = 3.")
	    ("n_queue_calls_eval", po::value<int>()->default_value(3),
	    	"In the Model Heuristic, number of calls in the queue that will be evaluated via optimization model. default = 3.")
	    ("n_scenarios", po::value<int>()->default_value(100),
	    	"Number of generated scenarios from distribution data.")
	    ("n_time_slots", po::value<int>()->default_value(12),
	    	"Number of time slots used in model. default = 12")
	    ("n_regions", po::value<int>()->default_value(76),
	    	"Number of regions used in model. default = 76")
	    ("n_hospitals", po::value<int>()->default_value(10),
	    	"Number of hospitals. default = 10, max number = 10")
	    ("n_bases", po::value<int>()->default_value(10),
	    	"Number of bases. default = 10, max number = 33")
	    ("n_cleaning_bases", po::value<int>()->default_value(3),
	    	"Number of cleaning bases. default = 3, max number = 10")
	    ("n_ambulances", po::value<int>()->default_value(10),
	    	"Number of ambulances. default = 10")
		("EPS", po::value<double>()->default_value(0.0001),
			"Tolerance")
		("min_preparedness", po::value<double>()->default_value(0.4),
			"Minimum preparedness level allowed for any region.")
	    ("osrm_map_path", po::value<std::string>()->default_value(
	    	"utils/osrm/RiodeJaneiro.osrm"),
	    	"Path for the .osrm file used for Open Street Map time and distance library.")
	    ("extended_model", po::bool_switch()->default_value(false),
	    	"In the Model, uses all sets of variables (don't assume near hospital/base). default = false.")
	    ("h_random_hospital", po::bool_switch()->default_value(false),
	    	"In Heuristic, calls go to a random hospital. default = false")
	;

	po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description config_file_options;
	config_file_options.add(config);

	po::options_description visible("Allowed Options");
	visible.add(generic).add(config);

	po::variables_map vm;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).run(), vm);
    notify(vm);

	std::ifstream ifs(config_file.c_str());
    if (vm.count("file") && !ifs){
        std::cout << "Unable to open file: " << config_file << "\n";
        exit(1);
    }
    else if(vm.count("file") && ifs){
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
    }

    if (vm.count("help")) {
        std::cout << visible << "\n";
    }
    if (vm.count("version")) {
    	//TODO: ler versão do SVN?
        version();
    }

    return Param(vm);
}



double f(int t, int c, int a, int b_h1, int l, int h, Data& data,
	int factor){
	// std::cout << t << " " << c << " " << a << " " << b_h1 << " ";
	// std::cout << l << " " << h << ", factor " << factor  << "\n";
	// return data.tao[t][c][a][b_h1][l][h]*factor*data.ins.penalty_matrix[a][c];
	double response_time = data.tao[t][c][a][b_h1][l][h] - 1200;
	return response_time*factor*data.ins.penalties[c] + 3600*data.ins.penalty_matrix[a][c];
}


double f(int t, int c, int a, int l1, int b, int l, int h,
	Data& data, int factor){
	// return data.tao[t][c][a][l1][l][h]*factor*data.ins.penalty_matrix[a][c];
	double response_time = data.tao[t][c][a][l1][l][h] - 1200;
	return response_time * factor * data.ins.penalties[c] + 3600*data.ins.penalty_matrix[a][c];
}

/*Cost functions of bases and queue.*/
double g_tcl(int t, int c, int l, double factor){
	return 1*factor;
}

double g_tab(int t, int a, int b, double factor){
	return 1*factor;
}

/*b must be the real location, not its index. You should pass bases[b] instead
of b*/
double g_talb(int t, int a, int l1, int b, Data& data, double factor){
	return 1*factor;
}


// double travel_time = osrm_helper.get_duration(std::make_pair(-22.900434, -43.277936),
// 	std::make_pair(-22.859857, -43.307751));
// std::cout << "Travel Time OSRM: " << travel_time << "\n";


	// std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);


	// GRBEnv env;
	// GeneratorRegressor gen(env, calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_old = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T regressor old(s): {}\n", run_time_quick);

	// QuickGeneratorRegressor qr_gen(env,calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// qr_gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_quick = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T regressor quick(s): {}\n", run_time_quick);

	// Generator gen(calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_old = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T No regressor (s): {}\n", run_time_old);


	// QuickGenerator q_gen(calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// q_gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_quick = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T no regressor quick(s): {}\n", run_time_quick);
	// fmt::print("T old(s): {}, T new(s): {}\n", run_time_old, run_time_quick);




	
	// GRBEnv env;
	// GeneratorRegressor gen(env, calls_path,neighbors_path, info_path);
	// gen.test();

	// fmt::print("Finished instance {}\n",g_params.generator_folder);