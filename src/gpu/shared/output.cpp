#include "output.hpp"

#include <stdexcept>

#include <cstdio>
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


namespace output
{
	Output load( char const* ) //XXX-stub; 
	{
		return {};
	}
	

#	if 0
	Output load( char const* aPath ) try
	{
		// Load main "output" XML
		tinyxml2::XMLDocument output;
		{
			auto xerr = output.LoadFile( aPath );

			if( tinyxml2::XML_SUCCESS != xerr )
			{
				throw error::XMLGeneric( "TinyXML: LoadFile: %s", tinyxml2::XMLDocument::ErrorIDToName(xerr) );
			}
		}

		auto const oroot = output.FirstChildElement( "output" );
		if( !oroot )
			throw error::XMLMissing( "No <output> root element" );

		// load data
		Output ret;

		try
		{
			ret.model = from_string(child_text_(oroot, "model" ), Identity<EModel>{});
		}
		catch( error::XMLLoadError const& )
		{
			bool ok = false;
			try
			{
				ret.model = from_string(child_text_(oroot, "distribution_class"), Identity<EModel>{});
				ok = true;

				std::fprintf( stderr, "WARNING: using deprecated \"distribution_class\"\n" );
			}
			catch( error::XMLLoadError const& )
			{}

			if( !ok ) throw;
		}
		
		ret.componentCount = child_int64_( oroot, "number_of_components" );
		ret.abcSampleCount = child_int64_( oroot, "number_of_abc_samples" );

		parse_list_( oroot, "m",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.m.push_back( aX ); }
		);
		parse_list_( oroot, "s",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.s.push_back( aX ); }
		);
		parse_list_( oroot, "c",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.c.push_back( aX ); }
		);
		parse_list_( oroot, "az",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.az.push_back( aX ); }
		);
		parse_list_( oroot, "dist",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.dist.push_back( aX ); }
		);
		parse_list_( oroot, "w",
			[] (char const* aX, char** aY) { return std::strtod( aX, aY ); },
			[&ret] (double aX) { ret.w.push_back( aX ); }
		);

		ret.epsilon = child_double_( oroot, "epsilon" );

		// optionally, load the meta data
		if( auto const source = oroot->FirstChildElement( "source" ) )
		{
			char const* appName = source->Attribute( "appName" );
			char const* date = source->Attribute( "date" );
			char const* machine = source->Attribute( "machine" );

			ret.producerName = appName ? appName : "<not specified>";
			ret.date = date ? date : "<not specified>";
			ret.machine = machine ? machine : "<not specified>";
		}
		else
		{
			ret.producerName = ret.date = ret.machine = "<unknown>";
		}

		return ret;
	}
	catch( ... )
	{
		std::throw_with_nested( error::XMLLoadError( "Error loading '%s'", aPath ) );
	}
#	endif
}


namespace output
{
	void write( char const* aPath, Output const& aOutput )
	{
		if( auto fof = std::fopen( aPath, "wb" ) )
		{
			std::fprintf( fof, "<output>\n" );

			// write simulation params
			std::fprintf( fof, "\t<model>%s</model>\n", to_string(aOutput.model).c_str() );
			std::fprintf( fof, "\t<weighting_scheme>%s</weighting_scheme>\n", to_string(aOutput.weightingScheme).c_str() );
			std::fprintf( fof, "\t<number_of_abc_samples>%zu</number_of_abc_samples>\n", aOutput.abcSampleCount );
			std::fprintf( fof, "\t<number_of_components>%zu</number_of_components>\n", aOutput.componentCount );
			std::fprintf( fof, "\t<number_of_z_components>%zu</number_of_z_components>\n", aOutput.zComponentCount );
			std::fprintf( fof, "\n" );

			// write meta data
			std::fprintf( fof, "\t<meta>\n" );
			for( auto it = aOutput.meta.begin(); it != aOutput.meta.end(); ++it )
			{
				std::fprintf( fof, "\t\t<%s>%s</%s>\n", it->first.c_str(), it->second.c_str(), it->first.c_str() );
			}
			std::fprintf( fof, "\t</meta>\n" );
			std::fprintf( fof, "\n" );

			// write results
			std::fprintf( fof, "\t<converged>%s</converged>\n", aOutput.converged ? "true" : "false" );
			std::fprintf( fof, "\t<epsilon>%.18g</epsilon>\n", aOutput.epsilon );
			std::fprintf( fof, "\t<m>" );
			{
				auto i = aOutput.m.begin();
				for( auto const  e = aOutput.m.end()-1; i != e; ++i )
					std::fprintf( fof, "%.18g, ", *i );
				std::fprintf( fof, "%.18g", *i );
			}
			std::fprintf( fof, "\t</m>\n" );

			std::fprintf( fof, "\t<c>" );
			{
				auto i = aOutput.c.begin();
				for( auto const e = aOutput.c.end()-1; i != e; ++i )
					std::fprintf( fof, "%.18g, ", *i );
				std::fprintf( fof, "%.18g", *i );
			}
			std::fprintf( fof, "\t</c>\n" );

			std::fprintf( fof, "\t<az>" );
			{
				auto i = aOutput.az.begin();
				for( auto const e = aOutput.az.end()-1; i != e; ++i )
					std::fprintf( fof, "%.18g, ", *i );
				std::fprintf( fof, "%.18g", *i );
			}
			std::fprintf( fof, "\t</az>\n" );

			std::fprintf( fof, "\t<dist>" );
			{
				auto i = aOutput.dist.begin();
				for( auto const e = aOutput.dist.end()-1; i != e; ++i )
					std::fprintf( fof, "%.18g, ", *i );
				std::fprintf( fof, "%.18g", *i );
			}
			std::fprintf( fof, "\t</dist>\n" );

			std::fprintf( fof, "\t<w>" );
			{
				auto i = aOutput.w.begin();
				for( auto const e = aOutput.w.end()-1; i != e; ++i )
					std::fprintf( fof, "%.18g, ", *i );
				std::fprintf( fof, "%.18g", *i );
			}
			std::fprintf( fof, "\t</w>\n" );

			std::fprintf( fof, "</output>\n" );
			std::fclose( fof );
		}
	}
}
