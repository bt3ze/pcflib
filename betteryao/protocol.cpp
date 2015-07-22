#include <iostream>
#include <log4cxx/propertyconfigurator.h>

//#include "Yao.h"
//#include "BetterYao4.h"
#include "BetterYao5.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	log4cxx::PropertyConfigurator::configure("log4cxx.conf");

	if (argc < 8)
	{
		std::cout << "Usage:" << std::endl
			<< "\tbetteryao [secu_param] [stat_param] [pcf_file] [input_file] [ip_server] [port_base] [mode]" << std::endl
			<< std::endl
			<< "[secu_param]: multiple of 8 but 128 at most" << std::endl
			<< "[stat_param]: multiple of the cluster size" << std::endl
			<< "  [pcf_file]: from the PCF compiler" << std::endl
			<< " [ip_server]: the IP (not domain name) of the IP exchanger" << std::endl
			<< " [port_base]: for the \"IP address in use\" hassles" << std::endl
			<< "      [mode]: 0=>honest-but-curious, 1=>malicious" << std::endl
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	EnvParams params;

	YaoBase *sys = 0;

	params.secu_param   = atoi(argv[1]);
	params.stat_param   = atoi(argv[2]);

	//params.circuit_file = argv[3];
	params.pcf_file     = argv[3];
	params.private_file = argv[4];
	params.ipserve_addr = argv[5];

	params.port_base    = atoi(argv[6]);

        // The special file that holds the inputs
        params.input_file = argv[4];

	switch(atoi(argv[7]))
	{
        case 0:
          // for now, also building BY5, but it will be compiled
          // with a different BY5 file (the hbc one)
          sys = new BetterYao5(params);
          break;
	case 1:
          sys = new BetterYao5(params);
	  break;

	default:
          std::cout << "Argument Failure" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
          break;
        }

	sys->start();
	delete sys; // delete MPI objects before MPI_Finalize()

	MPI_Finalize();

	return 0;
}
