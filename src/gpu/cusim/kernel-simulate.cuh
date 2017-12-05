#ifndef KERNEL_SIMULATE_CUH_AD83D236_3F02_4617_ACE4_A8BAEAB08009
#define KERNEL_SIMULATE_CUH_AD83D236_3F02_4617_ACE4_A8BAEAB08009

#include "dsinline.cuh"
#include "particle.cuh"

namespace cusim
{
	/** Simulation setup
	 *
	 * Configuration and values that are shared across multiple simulations.
	 */
	template<
		typename tCount,
		typename tScalar,
		EModel tModel,
		class tComponentCount
	>
	struct SimulateSetup
	{
		using Count = tCount;
		using Scalar = tScalar;

		using CompCount = tComponentCount;
		using Particle = Particle<Scalar,CompCount>;

		static constexpr EModel kModel = tModel;

		CompCount compCount;

		Count jobCount;
		Count kmin;

		Scalar fourDeltaT;

		Scalar halfLx, halfLy, halfLz;
		Scalar halfAx, halfAy;
	};

	/** Simulation run data
	 *
	 * Values valid for a single simulation run.
	 *
	 * \note The template parameters have to match those of `SimulateSetup`.
	 */
	template<
		typename tCount,
		typename tScalar,
		class tZCount,
		class tComponentCount
	>
	struct SimulateRun
	{
		using Count = tCount;
		using Scalar = tScalar;

		Count const* frames; // TODO: move to system setup?
		Count const* particles;

		DSInline<Scalar,tZCount> halfAz;
		DSInline<Scalar,tComponentCount> preCompProb;
		DSInline<Scalar,tComponentCount> randWalkStddev;
	};


	template< 
		class tSimulateSetup,
		class tSimulateRun,
		class tOutput,
		class tRandData
	>
	__global__
	void K_simulate_system( tSimulateSetup, tSimulateRun, tOutput, tRandData );
}

#include "kernel-simulate.inl"
#endif // KERNEL_SIMULATE_CUH_AD83D236_3F02_4617_ACE4_A8BAEAB08009
