/*-******************************************************* -- HEADER -{{{1- */
/*-	simulation_impl.cuh
 *
 * Provides SimulationT<>, a template with god-object-like ambitions. It
 * implements running a complete estimation on selected GPUs.
 *
 * TODO:
 *   - rename CountType to DeviceCount; probably don't need a HostCount...
 *   - normal distribution via sim_arg ?
 */
/*-***************************************************************** -}}}1- */

#ifndef SIMULATION_IMPL_CUH_FC0DF63C_0AB8_4BBA_A6C8_60BDAC5F91DE
#define SIMULATION_IMPL_CUH_FC0DF63C_0AB8_4BBA_A6C8_60BDAC5F91DE

#include <tuple>
#include <random>
#include <type_traits>

#include <cstddef>
#include <cstdint>

#include "support/named.hpp"
#include "support/queue.hpp"
#include "support/matrix.hpp"
#include "support/dsvalue.hpp"
#include "support/dsvector.hpp"

#include "shared/input.hpp"
#include "shared/common.hpp"

#include "cusim/particle.cuh"
#include "cusim/dsinline.cuh"
#include "cusim/histogram2d.cuh"

#include "curng/uniform.cuh"
#include "curng/engine-lcg.cuh"
#include "curng/normal-boxmuller.cuh"

#include "pool.cuh"
#include "random.cuh"
#include "simulation.hpp"

#define SIM_KERNEL_TIMINGS 1

#if SIM_KERNEL_TIMINGS
#	include "evpool.cuh"
#endif // ~ SIM_KERNEL_TIMINGS

namespace sim_arg
{
	/** Aspect: Simulation model
	 *
	 * Statically determined simulation model.
	 */
	template< EModel tModel >
	using Model = named::ValueArgument< struct Model_, EModel, EModel::discreteFixedZ, tModel >;

	/** Aspect: host scalar type
	 *
	 * Scalar type used for floating-point computations on the host. May be 
	 * either float or double (and possibly long double).
	 */
	template< typename tType > 
	using HostScalar = named::TypeArgument< struct HostScalar_, float, tType >;
	/** Aspect: device scalar type
	 *
	 * Scalar type used for floating-point computations on the device. May be 
	 * either float or double (and possibly long double).
	 */
	template< typename tType > 
	using DeviceScalar = named::TypeArgument< struct DeviceScalar_, float, tType >;

	
	/** Aspect: count type
	 *
	 * Type used for counting the number of trials, and in the histogram. The
	 * values are integral and larger than or equal to zero.
	 */
	template< typename tType >
	using CountType = named::TypeArgument< struct CountType_, std::uint32_t, tType >;

	/** Aspect: host random engine type
     *
	 * Host random number engine.
	 */
	template< class tRng >
	using HostRng = named::TypeArgument< struct HostRng_, std::mt19937, tRng >;
	/** Aspect: device random engine type
	 *
	 * Device random number engine
	 */
	template< class tRng >
	using DeviceRng = named::TypeArgument< struct DevRng_, curng::EngineLCG48_32, tRng >;

	/** Aspect: component count
	 *
	 * Component count. By default uses a dynamic value, derived from the input
	 * parameters. Can optionally be fixed to a static value with `StaticValue`.
	 *
	 * Static values one and two are further optimized.
	 */
	template< class tCCount >
	using ComponentCount = named::TypeArgument< struct CompCount_, DynamicValue<unsigned>, tCCount >;

	/** Aspect: Matrix layout
	 *
	 * Internal matrix layout (row or column major)
	 */
	template< class tMatLayout >
	using MatrixLayout = named::TypeArgument< struct MatLayout_, aspect::MatrixRowMajor, tMatLayout >;
};



template< class... tArgs >
class SimulationT final : public Simulation
{
	public:
		static constexpr EModel kModel = named::get_value<EModel, sim_arg::Model, tArgs...>::value;

		using HScalar = named::get_type_t<sim_arg::HostScalar, tArgs...>;
		using DScalar = named::get_type_t<sim_arg::DeviceScalar, tArgs...>;

		using Count = named::get_type_t<sim_arg::CountType, tArgs...>;

		using DevRng = named::get_type_t<sim_arg::DeviceRng, tArgs...>;
		using HostRng = named::get_type_t<sim_arg::HostRng, tArgs...>;

		using DeviceRandom = Random<
			DevRng,
			DScalar,
			curng::NormalBoxMuller,
			curng::UniformReal
		>;

		struct SystemSetup;
		struct SystemRunData;

	public:
		SimulationT( input::Parameters const&, HostRng&, SimulationConfig const& );
		~SimulationT();

		SimulationT( SimulationT const& ) = delete;
		SimulationT& operator= (SimulationT const) = delete;

	public:
		void run( SimHostRNG& ) override;
		output::Output output() override;

	private:
		using AbcCount_ = DynamicValue<unsigned>;
		using CCount_ = named::get_type_t<sim_arg::ComponentCount, tArgs...>;
		using ZCount_ = typename std::conditional<
			kModel == EModel::discreteFixedZ,
			StaticValue<unsigned,1>,
			CCount_
		>::type;

		using MatLayout_ = named::get_type_t<sim_arg::MatrixLayout, tArgs...>;

		template< typename tType > using CVec_ = DSVector<tType, CCount_>;
		template< typename tType > using ZVec_ = DSVector<tType, ZCount_>;
		template< typename tType > using AVec_ = DSVector<tType, AbcCount_>;
		
		template< typename tType > using CAMat_ = Matrix<
			tType,
			MatLayout_,
			CCount_, AbcCount_
		>;
		template< typename tType > using ZAMat_ = Matrix<
			tType,
			MatLayout_,
			ZCount_, AbcCount_
		>;

		struct Sample_
		{
			SimulationT* parent;

			CVec_<HScalar> mBis, cBis, powCBis;
			ZVec_<HScalar> azBis;
			HScalar distBis;

			CVec_<std::size_t> indices;
			CVec_<HScalar> mTmp, cTmp;
			ZVec_<HScalar> azTmp;

			DScalar* halfAz;
			DScalar* preCompProb;
			DScalar* randWalkStddev;

			std::vector<Count> particleCounts;

			int device; //TODO: split this off elsewhere???

			SystemRunData sampleRunData;
			DScalar* devResultDistance;
			DScalar* hostResultDistance;

#			if SIM_KERNEL_TIMINGS
			cudaEvent_t simStart, simStop;
			cudaEvent_t distStart, distStop;
			cudaEvent_t totalStart, totalStop;
#			endif // ~ SIM_KERNEL_TIMINGS
		};

		struct Result_
		{
			Sample_* sample;
			cudaError_t error;
		};


		struct CudaDevGlobal_
		{
			int device;
			std::size_t queueCount;
			
			cusim::Histogram2D<Count> reference;

			Pool<Count> particleCountPool;
			MappedPool<DScalar> resultDistancePool;

			Pool<DScalar> halfAzPool; //XXX-NOTE: only needed if ZCount_ is dynamic
			Pool<DScalar> preCompProbPool; //XXX-NOTE: only needed if CCount_ is dynamic
			Pool<DScalar> randWalkStddevPool; //XXX-NOTE: only needed if CCount_ is dyn.

			dim3 blocks, threads;
			std::size_t totalThreads;

			Count const* devFrameCounts;

#			if SIM_KERNEL_TIMINGS
			EvPool timeEvents;
#			endif // ~ SIM_KERNEL_TIMINGS
		};
		struct CudaDevQueue_
		{
			std::size_t devidx;
			cudaStream_t stream;

			cusim::Histogram2D<Count> result;
			typename DeviceRandom::GlobalData randomState;
		};

		using WeightingFun_ = void (SimulationT::*)();

	public:
		/* CUDA/NVCC needs these to be public: "A type that is defined inside a
		 * class and has private or protected access ("...") cannot be used in
		 * the template argument type of a __global__ function template
		 * instantiation, unless the class is local to a __device__ or
		 * __global__ function."
		 */
		struct SystemSetup
		{
			using value_type = DScalar;
			using count_type = Count;

			using CompCount = CCount_;
			using Particle = cusim::Particle<value_type,CompCount>;

			static constexpr EModel kModel = SimulationT::kModel;

			CompCount compCount;

			count_type jobCount;
			count_type kmin;

			value_type deltaT; //TODO: 4*deltaT??

			value_type halfLx, halfLy, halfLz;

			value_type halfAx;
			value_type halfAy;
		};

		struct SystemRunData
		{
			using value_type = DScalar;
			using count_type = Count;
			
			count_type const* frames; //TODO: could be SystemSetup. :-/
			count_type const* particles;
			
			cusim::DSInline<value_type,ZCount_> halfAz;
			cusim::DSInline<value_type,CCount_> preCompProb;
			cusim::DSInline<value_type,CCount_> randWalkStddev;
		};

	private:
		void prepare_( input::Parameters const&, HostRng&, SimulationConfig const& );

		void compute_tau_();
		auto compute_initial_gamma_( HostRng& ) -> HScalar;

		void weighting_scheme_pmc_standard_();
		void weighting_scheme_inv_dist_sq_();

		void job_queue_( std::size_t, Sample_& );

		void sample_dev_init_( Sample_&, CudaDevGlobal_&, CudaDevQueue_& );
		void sample_dev_clean_( Sample_&, CudaDevGlobal_& );

		static void CUDART_CB cuda_stream_callback_(
			cudaStream_t,
			cudaError_t,
			void*
		);

	private:
		AbcCount_ mAbcCount;
		CCount_ mComponentCount;
		ZCount_ mZCount;

		SystemSetup mSystemSetup;

		std::size_t mKMax;
		std::size_t mDEBinCount;
		HScalar mDEBinWidth;
		HScalar mVolumeFactor;

		CAMat_<HScalar> mM, mC;
		ZAMat_<HScalar> mAz;
		AVec_<HScalar> mDist, mW, mPreW;
		CVec_<HScalar> mTauM, mTauC;
		ZVec_<HScalar> mTauAz;
		
		AVec_<HScalar> mWStar; //TODO: mWStar only if pmcStandard;

		WeightingFun_ mWeightingSchemeFn;

		HScalar mDeltaT;
		HScalar mGamma, mDeltaGamma;

		HScalar mAvgTrialCount;

		HScalar mLowerM, mUpperM;
		HScalar mLowerC, mUpperC;
		HScalar mLowerAz, mUpperAz;

		AVec_<Sample_> mSamples;
		Queue<Result_> mResults;

		std::vector<CudaDevGlobal_> mDevGlobal;
		std::vector<CudaDevQueue_> mDevQueues;

		std::size_t mNextQueue = 0;

		bool mConverged;
		std::size_t mMaxIter;

		int mVerbosity;

		std::size_t mInfoIterations, mInfoSimulations;

#		if SIM_KERNEL_TIMINGS
		float timeSimTotal, timeSimCount;
		float timeDistTotal, timeDistCount;
		float timeTotalTotal, timeTotalCount;
#		endif // ~ SIM_KERNEL_TIMINGS
};

namespace detail
{
	[[noreturn]] void throw_logic_error_( char const* );
	[[noreturn]] void throw_invalid_gpuspec_( char const* );
}

#include "detail/simulation_impl.inl"
#endif // SIMULATION_IMPL_CUH_FC0DF63C_0AB8_4BBA_A6C8_60BDAC5F91DE
