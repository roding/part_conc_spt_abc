#include "input.hpp"

#include <stdexcept>

#include <cassert>
#include <cstdlib>

#include <tinyxml2/tinyxml2.h>
#include <tinyformat/tinyformat.h>

#include "errors.hpp"

namespace
{
	tinyxml2::XMLElement const* child_checked_( tinyxml2::XMLElement const* aParent, char const* aChild )
	{
		if( auto const child = aParent->FirstChildElement( aChild ) )
			return child;

		throw error::XMLMissing( "No <%s> in <%s>", aChild, aParent->Name() );
	}
	
	char const* child_text_( tinyxml2::XMLElement const* aParent, char const* aChild )
	{
		return child_checked_( aParent, aChild )->GetText();
	}
	std::int64_t child_int64_( tinyxml2::XMLElement const* aParent, char const* aChild )
	{
		std::int64_t ret;
		auto xerr = child_checked_( aParent, aChild )->QueryInt64Text( &ret );
	
		if( tinyxml2::XML_SUCCESS != xerr )
		{
			throw error::XMLGeneric( "TinyXML: QueryInt64Text: %s for <%s> in <%s>", tinyxml2::XMLDocument::ErrorIDToName(xerr), aChild, aParent->Name() );
		}

		return ret;
	}
	double child_double_( tinyxml2::XMLElement const* aParent, char const* aChild )
	{
		double ret;
		auto xerr = child_checked_( aParent, aChild )->QueryDoubleText( &ret );
	
		if( tinyxml2::XML_SUCCESS != xerr )
		{
			throw error::XMLGeneric( "TinyXML: QueryDoubleText: %s for <%s> in <%s>", tinyxml2::XMLDocument::ErrorIDToName(xerr), aChild, aParent->Name() );
		}

		return ret;
	}
	bool child_bool_( tinyxml2::XMLElement const* aParent, char const* aChild )
	{
		bool ret;
		auto xerr = child_checked_( aParent, aChild )->QueryBoolText( &ret );
	
		if( tinyxml2::XML_SUCCESS != xerr )
		{
			throw error::XMLGeneric( "TinyXML: QueryBoolText: %s for <%s> in <%s>", tinyxml2::XMLDocument::ErrorIDToName(xerr), aChild, aParent->Name() );
		}

		return ret;
	}


	input::Bounds get_bounds_( tinyxml2::XMLElement const* aParent, char const* aPattern )
	{
		std::string const lo = tfm::format( aPattern, "lb" );
		auto const lb = child_double_( aParent, lo.c_str() );

		std::string const up = tfm::format( aPattern, "ub" );
		auto const ub = child_double_( aParent, up.c_str() );

		return { lb, ub };
	}

	template< typename tParse, typename tSink >
	void parse_list_( tinyxml2::XMLElement const* aParent, char const* aChild, tParse const& aParse, tSink const& aSink )
	{
		char const* ptr = child_text_( aParent, aChild );
		char* next;
		for( auto i = aParse( ptr, &next ); ptr != next; i = aParse( ptr, &next ) )
		{
			ptr = next;
			if( ',' == *ptr ) 
				++ptr;
				
			aSink( i );
		}

		assert( ptr == next );
		if( *ptr != '\0' )
		{
			throw error::XMLParsing( "Can't parse list <%s> in <%s>. Unparsable mayhem '%s' remains.", aChild, aParent->Name(), ptr );
		}
	}
}

namespace input
{

	Parameters load( char const* aPath, char const* aDataPath ) try
	{
		// Load main input XML
		tinyxml2::XMLDocument input;
		{
			auto xerr = input.LoadFile( aPath );

			if( tinyxml2::XML_SUCCESS != xerr )
			{
				throw error::XMLGeneric( "TinyXML: LoadFile: %s (input)", tinyxml2::XMLDocument::ErrorIDToName(xerr) );
			}
		}

		auto const iroot = input.FirstChildElement( "input" );
		if( !iroot )
		{
			throw error::XMLMissing( "No <input> root element (input)" );
		}

		// Load "data" XML
		/*char const* dataPath = aDataPath
			? aDataPath
			: child_text_( iroot, "data_file_path" )
		;*/
		if( !aDataPath )
			aDataPath = child_text_( iroot, "data_file_path" );

		tinyxml2::XMLDocument data;
		{
			auto xerr = data.LoadFile( aDataPath );

			if( tinyxml2::XML_SUCCESS != xerr )
			{
				throw error::XMLGeneric( "TinyXML: LoadFile: %s (data)", tinyxml2::XMLDocument::ErrorIDToName(xerr) );
			}
		}

		auto const droot = data.FirstChildElement( "data" );
		if( !droot )
			throw error::XMLMissing( "No <data> root element (data)" );
		
		// Extract settings
		Parameters ret;

		// input.xml
		char const* model = nullptr;
		try
		{
			model = child_text_( iroot, "model" );
		}
		catch( std::runtime_error const& )
		{
			try 
			{ 
				model = child_text_( iroot, "distribution_class" ); 
				std::fprintf( stderr, "WARNING: using deprecated \"distribution_class\" field\n" );
			}
			catch( std::runtime_error const& )
			{}

			if( !model )
				throw;
		}
		
		ret.model            = from_string( model, Identity<EModel>{} );

		char const* wscheme = child_text_( iroot, "weighting_scheme" );
		ret.weightingScheme  = from_string( wscheme, Identity<EWeightingScheme>{} );

		ret.componentCount   = child_int64_( iroot, "number_of_components" );
		ret.kmin             = child_int64_( iroot, "kmin" );
		ret.deBinCount       = child_int64_( iroot, "number_of_de_bins" );
		ret.abcSampleCount   = child_int64_( iroot, "number_of_abc_samples" );
		ret.avgTrialCount    = child_int64_( iroot, "ub_average_number_of_trials" );

		ret.Lx               = child_double_( iroot, "Lx" );
		ret.Ly               = child_double_( iroot, "Ly" );
		ret.Lz               = child_double_( iroot, "Lz" );
		ret.gamma            = child_double_( iroot, "gamma_initial" );
		ret.deltaGamma       = child_double_( iroot, "delta_gamma" );

		try
		{
			ret.adaptiveGamma = child_bool_( iroot, "gamma_adaptive" );
		}
		catch( std::exception const& )
		{
			std::fprintf( stderr, "Note: <gamma_adaptive> not set/invalid (old input?); defaulting to false.\n" );
			ret.adaptiveGamma = false;
		}

		// ret.de seems special: the input only defines the upper bound.
		ret.de.lower         = 0.f;
		ret.de.upper         = child_double_( iroot, "ub_de" );

		ret.m                = get_bounds_( iroot, "%s_m" );
		try
		{
			ret.s = get_bounds_( iroot, "%s_s" );
			std::fprintf( stderr, "Note: deprecated <*_s> tags found (old input?)\n" );
		}
		catch( error::XMLMissing const& )
		{
			//ignore
		}
		ret.c                = get_bounds_( iroot, "%s_c" );
		ret.az               = get_bounds_( iroot, "%s_az" );

		ret.outputFilePath   = child_text_( iroot, "output_file_path" );

		// data.xml
		ret.ax               = child_double_( droot, "ax" );
		ret.ay               = child_double_( droot, "ay" );
		ret.deltaT           = child_double_( droot, "deltat" );
		
		parse_list_( droot, "number_of_frames",
			[] (char const* aX, char** aY) { return std::strtoull( aX, aY, 0 ); },
			[&ret] (std::size_t aX) { ret.frameCounts.push_back( aX ); }
		);
		parse_list_( droot, "K",
			[] (char const* aX, char** aY) { return std::strtoull( aX, aY, 0 ); },
			[&ret] (std::size_t aX) { ret.Ks.push_back( aX ); }
		);
		parse_list_( droot, "DE",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.DEs.push_back( aX ); }
		);

		return ret;
	}
	catch( ... )
	{
		std::throw_with_nested( error::XMLLoadError( "Error loading input=%s and data=%s", aPath, aDataPath ) );
	}
}
