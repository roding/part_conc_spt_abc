#include <chrono>
#include <string>
#include <vector>
#include <typeinfo>
#include <stdexcept>

#include <cstdio>
#include <cassert>
#include <cstddef>
#include <cstdint>
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
		EAction_ action      = EAction_::compute;
		
		std::string input    = "./input.xml";  // main input "input.xml"
		std::string data;      // if set, overrides the data.xml path from the input
		std::string output;    // if set, overrides the output.xml path from the input

		ESimScalar scalar    = ESimScalar::floatType;

		int verbosity        = 0;
		
		std::string gpuSpec  = "1/30";

		std::size_t maxIter  = 0;

		bool fixedSeed       = false;
		std::uint32_t seedValue;

		std::string self;
	};

	void help_( FILE*, char const* );
	Config_ parse_command_line_( int, char*[] );

	input::Parameters override_params_( input::Parameters, Config_ const& );
}

int main( int aArgc, char* aArgv[] ) try
{
	using SysClock_ = std::chrono::system_clock;
	
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

			SimulationConfig scfg;
			scfg.scalarType  = cfg.scalar;
			scfg.maxIter     = cfg.maxIter;
			scfg.gpuSpec     = cfg.gpuSpec;
			scfg.verbosity   = cfg.verbosity;

			auto rng = [&] () {
				if( cfg.fixedSeed )
					return SimHostRNG( cfg.fixedSeed );

				return make_prng<SimHostRNG>();
			}();
			auto sim = make_simulation( rng, param, scfg );

			auto const start = SysClock_::now();
			sim->run( rng );
			auto const end = SysClock_::now();


			auto output = sim->output();
			{
				output.meta["runtime_ms"] = tfm::format( "%u", std::chrono::duration_cast<std::chrono::duration<std::uint64_t,std::milli>>(end-start).count() );

				// no put_time() on GCC. :-(
				auto startTime = SysClock_::to_time_t(start);
				auto endTime = SysClock_::to_time_t(end);
				
				char buff[256];
				std::strftime( buff, 255, "%F %T%z", std::localtime(&startTime) );
				output.meta["started_at"] = buff;

				std::strftime( buff, 255, "%F %T%z", std::localtime(&endTime) );
				output.meta["finished_at"] = buff;

#				ifdef NDEBUG
				output.meta["debug"] = "false";
#				else
				output.meta["debug"] = "true";
#				endif
			}
			
			output::write( param.outputFilePath.c_str(), output );
			break;
		}
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

  --gpus, -g <gpuspec>  : select GPUs and queues
  --scalar, -s {f,d}    : select scalar type (𝗳loat or 𝗱ouble)

  --max-steps, -S <N>   : abort after <N> steps
  --fixed-seed, -F <N>  : use a fixed seed <N> for random number generation

  --quiet, -q           : decrease verbosity
  --verbose, -v         : increase verbosity (repeat for further increases)
  
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
			else if( flag_( "--quiet", "-q" ) )
			{
				--cfg.verbosity;
			}
			else if( flag_( "--verbose", "-v" ) )
			{
				++cfg.verbosity;
			}
			else if( arg_( "--input", "-i", cfg.input ) )
			{}
			else if( arg_( "--data", "-d", cfg.data ) )
			{}
			else if( arg_( "--output", "-o", cfg.output ) )
			{}	
			else if( arg_( "--gpus", "-g", cfg.gpuSpec ) )
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
			else if( flag_( "--max-steps", "-S" ) )
			{
				if( arg+1 >= aArgc ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional integer argument", aArgv[arg] ) );

				char dummy;
				int iret = std::sscanf( aArgv[arg+1], "%zu%c", &cfg.maxIter, &dummy );

				if( 1 != iret ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional positive integer argument and not %s", aArgv[arg], aArgv[arg+1] ) );

				++arg;
			}
			else if( flag_( "--fixed-seed", "-F" ) )
			{
				if( arg+1 >= aArgc ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional integer argument", aArgv[arg] ) );

				char dummy;
				int iret = std::sscanf( aArgv[arg+1], "%u%c", &cfg.seedValue, &dummy );

				if( 1 != iret ) throw std::runtime_error( tfm::format( "Command line argument %s requires an additional integer argument and not %s", aArgv[arg], aArgv[arg+1] ) );

				cfg.fixedSeed = true;

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
