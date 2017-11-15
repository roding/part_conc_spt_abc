#ifndef PARTICLE_CUH_F641B46F_9145_43A1_BD4A_DFD20ADAE317
#define PARTICLE_CUH_F641B46F_9145_43A1_BD4A_DFD20ADAE317

#include <type_traits>

#include "../support/dsvalue.hpp"

namespace cusim
{
	template< typename tReal, class tCompCount >
	struct Particle
	{
		tReal x, y, z;
		unsigned index;

		static_assert( 
			std::is_same<typename tCompCount::type, unsigned>::value,
			"cusim: use 'unsigned' component counts"
		);
	};

	template< typename tReal >
	struct Particle< tReal, StaticValue<unsigned,1> >
	{
		tReal x, y, z;
	};
}

#endif // PARTICLE_CUH_F641B46F_9145_43A1_BD4A_DFD20ADAE317
