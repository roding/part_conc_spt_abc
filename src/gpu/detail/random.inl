#include "../support/utility.hpp"

#define RANDOM_TEMPL_ template< class tEng, typename tReal, class tNDist, class tUDist >
#define RANDOM_CLASS_ Random<tEng,tReal,tNDist,tUDist>

RANDOM_TEMPL_ __device__ inline
RANDOM_CLASS_::Instance::Instance( unsigned aTid, GlobalData& aData )
	: tEng::Engine( aTid, aData )
	, tNDist::template Distribution<tReal>( aTid, aData )
	, tUDist::template Distribution<tReal>( aTid, aData )
{}

RANDOM_TEMPL_ __device__ inline
tReal RANDOM_CLASS_::Instance::uniform01( unsigned aTid, GlobalData& aData )
{
	using Eng_ = typename tEng::Engine;
	using UDist_ = typename tUDist::template Distribution<tReal>;
	return static_cast<UDist_&>(*this)( static_cast<Eng_&>(*this), aTid, aData, aData );
}
RANDOM_TEMPL_ __device__ inline
tReal RANDOM_CLASS_::Instance::normal01( unsigned aTid, GlobalData& aData )
{
	using Eng_ = typename tEng::Engine;
	using NDist_ = typename tNDist::template Distribution<tReal>;
	return static_cast<NDist_&>(*this)( static_cast<Eng_&>(*this), aTid, aData, aData );
}

RANDOM_TEMPL_ __device__ inline
void RANDOM_CLASS_::Instance::store( unsigned aTid, GlobalData& aData )
{
	static_cast<typename tEng::Engine&>(*this).store( aTid, aData );
}

RANDOM_TEMPL_ template< class tHostRNG > __host__ inline
auto RANDOM_CLASS_::initialize( std::size_t aThreadCount, tHostRNG& aHostRNG ) -> GlobalData
{
	using EngData_ = typename tEng::GlobalData;
	using NData_ = typename tNDist::template GlobalData<tReal>;
	using UData_ = typename tUDist::template GlobalData<tReal>;
	
	GlobalData ret;
	static_cast<EngData_&>(ret) = tEng::initialize( aThreadCount, aHostRNG );
	static_cast<NData_&>(ret) = tNDist::initialize( Identity<tReal>{}, aThreadCount );
	static_cast<UData_&>(ret) = tUDist::initialize( Identity<tReal>{}, aThreadCount );
	return ret;
}
RANDOM_TEMPL_ __host__ inline
void RANDOM_CLASS_::cleanup( GlobalData& aData )
{
	tEng::cleanup( aData );
	tNDist::cleanup( aData );
	tUDist::cleanup( aData );
}

#undef RANDOM_CLASS_
#undef RANDOM_TEMPL_

