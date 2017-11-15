#include <string>
#include <vector>
#include <typeinfo>
#include <stdexcept>

#include <cstdio>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include <tinyformat/tinyformat.h>

#include "sysinfo.hpp"
#include "simulation.hpp"

#include "shared/input.hpp"

namespace
{
	enum class EAction_
	{
		compute,
		listDevices
	};
	
	struct Config_
	{
		EAction_ action     = EAction_::compute;
		
		std::string input   = "./input.xml";  // main input "input.xml"
		std::string data;     // if set, overrides the data.xml path from the input
		std::string output;   // if set, overrides the output.xml path from the input

		ESimScalar scalar   = ESimScalar::floatType;
		
		std::string gpuSpec;

		std::string self;
	};

	void help_( FILE*, char const* );
	Config_ parse_command_line_( int, char*[] );

	input::Parameters override_params_( input::Parameters, Config_ const& );
}

int main( int aArgc, char* aArgv[] ) try
{
	Config_ const cfg = parse_command_line_( aArgc, aArgv );

	switch( cfg.action )
	{
		case EAction_::listDevices:
		{
			list_devices();
			break;
		}

		case EAction_::compute:
		{
			input::Parameters const param = override_params_(
				input::load( 
					cfg.input.c_str(), 
					cfg.data.empty() ? nullptr : cfg.data.c_str() 
				),
				cfg
			);

			//TODO: select devices

			SimulationConfig scfg;
			scfg.scalarType = cfg.scalar;
			scfg.gpuSpec = cfg.gpuSpec;

			auto rng = make_prng<SimHostRNG>();
			auto sim = make_simulation( rng, param, scfg );

			sim->run( rng );
			sim->write_results( param );

			break;
		}

		default: assert(false);
	}

	return 0;
}
catch( std::exception const& eErr )
{
	fprintf( stderr, "Top-level exception: %s\n", typeid(eErr).name() );
	fprintf( stderr, "  - %s\n", eErr.what() );

	try
	{
		std::rethrow_if_nested(eErr);
	}
	catch( std::exception const& eNest )
	{
		fprintf( stderr, "  - nested exception %s\n", typeid(eNest).name() );
		fprintf( stderr, "    - %s\n", eNest.what() );
	}

	return 1;
}


namespace
{
	void help_( FILE* aStream, char const* aSelfName )
	{
		static char const* const kFlagText = R"(Flags:
  --input, -i <file>    : select input XML file
  --data, -d <file>     : override data XML file location(*)
  --output, -o <file>   : override output XML file destination(*)

  --scalar, -s {f,d}    : select scalar type (𝗳loat or 𝗱ouble)

  --gpus, -g <gpuspec>  : select GPUs and queues

  --help, -h            : print this help and exit
  --list-devices        : list compute devices and exit
)";

		static char const* const kNotesText = R"(
<gpuspec> is a comma separated list of GPU IDs (as shown by --list-devices).
Optionally, a queue count can be added to each GPU ID by appending a "/N".
For example "-g 1/5,2" would select GPUs 1 and 2, with 5 queues on GPU 1."

(*) If specified, the -d/-o options override the paths given in the input XML
    for the data and for the output XML, respectively.
)";

		
		std::fprintf( aStream, "Synopsis: %s [--flags...]\n", aSelfName );
		std::fprintf( aStream, "\n" );
		std::fprintf( aStream, "%s", kFlagText );
		std::fprintf( aStream, "%s", kNotesText );
	}
}
namespace
{
	Config_ parse_command_line_( int aArgc, char* aArgv[] )
	{
		Config_ cfg;
		cfg.self = aArgv[0];

		for( int arg = 1; arg < aArgc; ++arg )
		{
			auto flag_ = [&] (char const* aMain, char const* aAlt = nullptr ) {
				return 0 == std::strcmp( aMain, aArgv[arg] ) || (aAlt && 0 == std::strcmp( aAlt, aArgv[arg] ));
			};
			auto arg_ = [&] (char const* aMain, char const* aAlt, std::string& aOut) {
				if( flag_( aMain, aAlt ) )
				{
					if( arg+1 >= aArgc ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional argument", aArgv[arg] ) );

					aOut = aArgv[++arg];
					return true;
				}
				
				return false;
			};
			
			if( flag_( "--help", "-h" ) )
			{
				help_( stdout, aArgv[0] );
				std::exit( 0 );
			}
			else if( arg_( "--input", "-i", cfg.input ) )
			{}
			else if( arg_( "--data", "-d", cfg.data ) )
			{}
			else if( arg_( "--output", "-o", cfg.output ) )
			{}
			else if( flag_( "--scalar", "-s" ) )
			{
				if( arg+1 >= aArgc ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional argument f or d", aArgv[arg] ) );

				if( 0 == strcmp( "f", aArgv[arg+1] ) ) 
					cfg.scalar = ESimScalar::floatType;
				else if( 0 == strcmp( "d", aArgv[arg+1] ) ) 
					cfg.scalar = ESimScalar::doubleType;
				else throw std::runtime_error( tfm::format( "Scalar type must be either 'f' (float) or 'd' (double) and not '%s'", aArgv[arg+1] ) );

				++arg;
			}
			else if( flag_( "--gpus", "-g" ) )
			{
				if( arg+1 >= aArgc ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional argument", aArgv[arg] ) );

				cfg.gpuSpec = aArgv[arg+1];
				++arg;
			}
			else if( flag_( "--list-devices" ) ) cfg.action = EAction_::listDevices;
			else 
			{
				throw std::runtime_error( tfm::format( "Unknown command line argument '%s'", aArgv[arg] ) );
			}
		}

		return cfg;
	}
}

namespace
{
	input::Parameters override_params_( input::Parameters aParam, Config_ const& aCfg )
	{
		if( !aCfg.output.empty() )
			aParam.outputFilePath = aCfg.output;

		return aParam;
	}	
}
