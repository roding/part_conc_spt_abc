#include <cmath>
#include <cstdio>

#include "../gpu/shared/output.hpp"

namespace
{
	template< typename tIter >
	double avg_( tIter aBeg, tIter aEnd )
	{
		double count = 0, accum = 0;
		for( ; aBeg != aEnd; ++aBeg, ++count )
			accum += *aBeg;

		return accum / count;
	}
}

int main( int aArgc, char* aArgv[] ) try
{
	if( aArgc != 3 )
	{
		std::fprintf( stderr, "Error: need exactly two command line args.\n" );
		return 1;
	}

	auto a = output::load( aArgv[1] );
	auto b = output::load( aArgv[2] );

	std::printf( "A: by %s on %s (%s)\n", a.producerName.c_str(), a.machine.c_str(), a.date.c_str() );
	std::printf( "B: by %s on %s (%s)\n", b.producerName.c_str(), b.machine.c_str(), b.date.c_str() );

	auto am = avg_( a.m.begin(), a.m.end() );
	auto bm = avg_( b.m.begin(), b.m.end() );
	std::printf( "m    : %18g vs %18g (δ = %18g)\n", am, bm, std::abs(am-bm) );

	auto as = avg_( a.s.begin(), a.s.end() );
	auto bs = avg_( b.s.begin(), b.s.end() );
	std::printf( "s    : %18g vs %18g (δ = %18g)\n", as, bs, std::abs(as-bs) );

	auto ac = avg_( a.c.begin(), a.c.end() );
	auto bc = avg_( b.c.begin(), b.c.end() );
	std::printf( "c    : %18g vs %18g (δ = %18g)\n", ac, bc, std::abs(ac-bc) );
	
	auto aaz = avg_( a.az.begin(), a.az.end() );
	auto baz = avg_( b.az.begin(), b.az.end() );
	std::printf( "az   : %18g vs %18g (δ = %18g)\n", aaz, baz, std::abs(aaz-baz) );

	auto ad = avg_( a.dist.begin(), a.dist.end() );
	auto bd = avg_( b.dist.begin(), b.dist.end() );
	std::printf( "dist : %18g vs %18g (δ = %18g)\n", ad, bd, std::abs(ad-bd) );

	auto aw = avg_( a.w.begin(), a.w.end() );
	auto bw = avg_( b.w.begin(), b.w.end() );
	std::printf( "w    : %18g vs %18g (δ = %18g)\n", aw, bw, std::abs(aw-bw) );

	std::printf( "eps  : %18g vs %18g (δ = %18g)\n", a.epsilon, b.epsilon, std::abs(a.epsilon-b.epsilon) );

	return 0;
}
catch( std::runtime_error const& eErr )
{
	std::fprintf( stderr, "Top-level exception: %s\n", typeid(eErr).name() );
	std::fprintf( stderr, "  - %s\n", eErr.what() );

	try
	{
		std::rethrow_if_nested(eErr);
	}
	catch( std::runtime_error const& eNest )
	{
		std::fprintf( stderr, "  - nested exception %s\n", typeid(eNest).name() );
		std::fprintf( stderr, "    - %s\n", eNest.what() );
	}

	return 1;
}
