/*-******************************************************* -- HEADER -{{{1- */
/*-	random.cuh
 *
 * Ties together the random engine and distributions into a single object that
 * can (more) easily be passed to a CUDA kernel.
 */
/*-***************************************************************** -}}}1- */

#ifndef RANDOM_CUH_735F0FB5_46FB_42EB_8BC5_2AC099DBFD4B
#define RANDOM_CUH_735F0FB5_46FB_42EB_8BC5_2AC099DBFD4B

/** Helper: Random number generation
 *
 * Wraps a random engine \a tEngine, a normal distribution \a tNormalDist and
 * an uniform distribution \a tUniformDist into a single class. Generates
 * floating point numbers of type \a tReal.
 *
 * Usage:
   \code
   // host:
   using Random_ = Random< EngineLCG48_32, float, NormalBoxMuller, UniformReal >;
  
   auto hostRNG = ...;
   auto randomData = Random_::initialize( maxKernelThreadCount, hostRNG );
  
   ...
   kernel<<<...>>>( randomData );
   ...
   
   Random_::cleanup( randomData );
  
   // device:
   template< class tRandom, class tRandFactory = typename tRandom::Factory >
   __global__ void kernel( tRandom& aRandData, ... )
   {
       auto const tid = blockIdx.x*blockDim.x + threadIdx.x; // global thread index
       auto random = typename tRandFactory::Instance( tid, aRandData );
  
       ...
       float u01 = random.uniform01( tid, aRandData );
       float n01 = random.normal01( tid, aRandData );
       ...
  
       random.store( tid, aRandData );
   }
   \endcode
 *
 * foo
 */
template< class tEngine, typename tReal, class tNormalDist, class tUniformDist >
struct Random
{
	using result_type = tReal;
	
	using Engine = tEngine;
	using NormalDistribution = tNormalDist;
	using UniformDistribution = tUniformDist;

	/** Global state
	 *
	 * `GlobalData` holds state necessary to create per-thread random number
	 * generator instances. It may be passed to a CUDA kernel as an argument.
	 * `GlobalData` is created using the `initialize()` method.
	 */
	struct GlobalData
		: tEngine::GlobalData
		, tNormalDist::template GlobalData<tReal>
		, tUniformDist::template GlobalData<tReal>
	{
		using Factory = Random;
		using result_type = tReal;
	};
	
	/** Per-thread instance
	 *
	 * Per-thread random generator instance. Each thread owns its own copy of
	 * `Instance`, where necessary run-time state is stored.
	 *
	 * The per-thread state is loaded from the `GlobalData` via the `Instance`
	 * constructor. After use, the global state may be updated the `store()`
	 * method.
	 */
	struct Instance
		: tEngine::Engine
		, tNormalDist::template Distribution<tReal>
		, tUniformDist::template Distribution<tReal>
	{
		/** Constructor: initialize per-thread instance
		 *
		 * Initialize per-thread random generator instance. \a aTid must be a
		 * globally unique per-thread index, that does not exceed the maximal
		 * number of threads specified on initialization. \a aRandData contains
		 * the global state previously created with `initialize()` on the host.
		 *
		 * \note Both \a aTid and \a aRandData need to be passed to all methods
		 * of `Instance`. They are <b>not</b> cached by the object. 
		 */
		__device__
		Instance( unsigned aTid, GlobalData& aRandData );

		/** Return uniformly distributed random value ∈ [0,1)
		 */
		__device__
		tReal uniform01( unsigned, GlobalData& );
		/** Return normally distributed random value
		 *
		 * Return normally distributed random value with mean zero and unit
		 * deviation.
		 */
		__device__
		tReal normal01( unsigned, GlobalData& );

		
		/** Update the global state
		 *
		 * Update the global state from this instance's local state. This
		 * ensures that subsequent uses (in e.g., later kernel invocations)
		 * continue to generate new random numbers, rather than repeating the
		 * sequence observed by the current kernel/thread.
		 */
		__device__
		void store( unsigned, GlobalData& );
	};

	/** Initialize global data
	 *
	 * Allocate and initialize global state. Seeds the states of the CUDA 
	 * random engines using the host random number generator \a tHostRNG.
	 *
	 * State for up to \a std::size_t threads is allocated. The thread index
	 * (passed to `Instance`'s methods) must not exceed this number.
	 */
	template< class tHostRNG > __host__
	static GlobalData initialize( std::size_t, tHostRNG& );

	/** Clean up global data
	 *
	 * Frees any global state associated with the \a GlobalData.
	 */
	__host__
	static void cleanup( GlobalData& );
};

#include "detail/random.inl"
#endif // RANDOM_CUH_735F0FB5_46FB_42EB_8BC5_2AC099DBFD4B
