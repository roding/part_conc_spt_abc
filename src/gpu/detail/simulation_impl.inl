#include <chrono>
#include <numeric>
#include <utility>
#include <iterator>
#include <algorithm>

#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>

#include <unistd.h>

#include <tinyformat/tinyformat.h>

#include "../support/cuda_error.hpp"

#include "../cusim/kernel-distance.cuh"
#include "../cusim/kernel-simulate.cuh"

namespace detail
{
	// "uncorrected" 1/N variance
	template< typename tForwardIt > inline
	auto var_( tForwardIt aBeg, tForwardIt aEnd )
		-> typename std::iterator_traits<tForwardIt>::value_type
	{
		using Ty_ = typename std::iterator_traits<tForwardIt>::value_type;

		Ty_ accum = Ty_(0), count = Ty_(0);
		for( tForwardIt it = aBeg; it != aEnd; ++it, ++count )
			accum += *it;

		Ty_ const mean = accum / count;
		Ty_ var = Ty_(0);
		for( tForwardIt it = aBeg; it != aEnd; ++it )
		{
			auto const x = *it - mean;
			var += x*x;
		}

		return var / count;
	}

	// mean_()
	template< typename tScalar, typename tInputIt > inline
	auto mean_( tInputIt aBeg, tInputIt aEnd )
		-> decltype(std::declval<tScalar>()+*std::declval<tInputIt>())
	{
		using Ty_ = decltype(std::declval<tScalar>()+*std::declval<tInputIt>());

		Ty_ accum = Ty_(0), count = Ty_(0);
		for( ; aBeg != aEnd; ++aBeg, ++count )
			accum += *aBeg;

		return accum / count;
	}

	// normpdf_()
	template< typename tScalar > inline
	tScalar normpdf_( tScalar aX, tScalar aMu, tScalar aSigma )
	{
		auto const xm = aX - aMu;
		return std::exp( -tScalar(0.5)*xm*xm/(aSigma*aSigma) ) / (std::sqrt(tScalar(2) * tScalar(M_PI)) * aSigma);
	}

	// wrand_idx_(): choose random element by weight and return its index
	template< class tRng, typename tRandomIt > inline
	std::size_t wrand_idx_( tRng& aRng, tRandomIt aBeg, tRandomIt aEnd )
	{
		using Scalar_ = typename std::iterator_traits<tRandomIt>::value_type;
		Scalar_ const rr = std::uniform_real_distribution<Scalar_>{0,1}(aRng);
		auto const it = std::lower_bound( aBeg, aEnd, rr );
		return std::min( std::distance(aBeg,aEnd)-1, std::distance(aBeg,it) );
	}

	// displace_()
	template< class tRng, typename tScalar > inline
	tScalar displace_( tRng& aRng, tScalar aCurrent, tScalar aTau, tScalar aLower, tScalar aUpper )
	{
		std::normal_distribution<tScalar> ndist( 0.0, aTau );
		tScalar const delta = ndist( aRng );
		tScalar bis = aCurrent + delta;

		while( bis < aLower || bis > aUpper )
		{
			if( bis < aLower ) bis += tScalar(2) * (aLower-bis);
			if( bis > aUpper ) bis -= tScalar(2) * (bis-aUpper);
		}

		return bis;
	}

	// mat2vec_()
	template< typename tOut, class tMatrix > inline
	std::vector<tOut> mat2vec_( tMatrix const& aMat )
	{
		std::vector<tOut> ret;
		ret.reserve( aMat.n() * aMat.m() );

		for( std::size_t i = 0; i < aMat.m(); ++i )
		{
			auto const c = aMat.col( i );
			for( auto x : c )
				ret.push_back( x );
		}
		
		return ret;
	}
}


template< class... tArgs > inline
SimulationT<tArgs...>::SimulationT( input::Parameters const& aPar, HostRng& aRng, SimulationConfig const& aCfg )
	: mAbcCount( aPar.abcSampleCount )
	, mComponentCount( aPar.componentCount )
	, mZCount( EModel::discreteFixedZ == aPar.model ? 1 : aPar.componentCount )
	, mSystemSetup{ mComponentCount }
	, mM( mComponentCount, mAbcCount )
	, mC( mComponentCount, mAbcCount )
	, mAz( mZCount, mAbcCount )
	, mConverged(false)
{
	resize( mDist, mAbcCount, Scalar(0) );
	resize( mW, mAbcCount );
	resize( mPreW, mAbcCount );
	resize( mWStar, mAbcCount );

	resize( mSamples, mAbcCount );

	resize( mTauM, mComponentCount );
	resize( mTauC, mComponentCount );
	resize( mTauAz, mZCount );

	//mActiveSamples.reserve( mAbcCount );

	prepare_( aPar, aRng, aCfg );
}
template< class... tArgs > inline
SimulationT<tArgs...>::~SimulationT()
{
	for( auto& sample : mSamples )
	{
		release_host_ptr( sample.sampleRunData.halfAz, sample.halfAz );
		release_host_ptr( sample.sampleRunData.preCompProb, sample.preCompProb );
		release_host_ptr( sample.sampleRunData.randWalkStddev, sample.randWalkStddev );
	}

	for( auto& queue : mDevQueues )
	{
		CUDA_CHECKED cudaSetDevice( mDevGlobal[queue.devidx].device );

		queue.result.free_cuda_resources();
		DeviceRandom::cleanup( queue.randomState );

		CUDA_CHECKED cudaStreamDestroy( queue.stream );
	}
	for( auto& dev : mDevGlobal )
	{
		CUDA_CHECKED cudaSetDevice( dev.device );

		dev.reference.free_cuda_resources();
		dev.particleCountPool.free_cuda_resources();
		dev.resultDistancePool.free_cuda_resources();
		dev.halfAzPool.free_cuda_resources();
		dev.preCompProbPool.free_cuda_resources();
		dev.randWalkStddevPool.free_cuda_resources();

		cudaFree( const_cast<Count*>(dev.devFrameCounts) );

#		if SIM_KERNEL_TIMINGS
		dev.timeEvents.free_cuda_resources();
#		endif // ~ SIM_KERNEL_TIMINGS
	}
}

template< typename... tArgs > inline
void SimulationT<tArgs...>::run( SimHostRNG& aRng )
{
	using Clock_ = std::chrono::high_resolution_clock;

	mConverged = false;
	mInfoIterations = mInfoSimulations = 0;

	
	std::vector<std::size_t> waitingSamples( mAbcCount );

	AVec_<Count> trials;
	resize( trials, mAbcCount );

	Scalar infoLambdaAcc, infoLambdaCount;
	Clock_::time_point infoIterStart, infoIterEnd;

	if( mGamma < Scalar(0) )
	{
		if( mVerbosity >= 0 ) std::printf( "Compute intial γ...\n" );
		mGamma = compute_initial_gamma_( aRng );
		if( mVerbosity >= 0 ) std::printf( "  initial γ = %g\n", mGamma );
	}

	if( mVerbosity >= 0 ) std::printf( "Begin main iterations...\n" );
	while( !mConverged )
	{
		++mInfoIterations;

		if( mMaxIter != 0 && mInfoIterations >= mMaxIter )
		{
			std::printf( "NOTE: --max-steps reached. stopping.\n" );
			break;
		}

		mGamma -= mDeltaGamma;

		std::fill_n( trials.begin(), mAbcCount, Count(0) );
		std::partial_sum( mW.begin(), mW.end(), mPreW.begin() );

		waitingSamples.resize( mAbcCount );
		std::iota( waitingSamples.begin(), waitingSamples.end(), 0 );

		for( auto& sample : mSamples )
			sample.distBis  = std::numeric_limits<Scalar>::infinity();

		
		infoIterStart = Clock_::now();
		infoLambdaAcc = infoLambdaCount = Scalar(0);

#		if SIM_KERNEL_TIMINGS
		timeSimCount = timeDistCount = timeTotalCount = 0.0f;
		timeSimTotal = timeDistTotal = timeTotalTotal = 0.0f;
#		endif // ~ SIM_KERNEL_COUNT

		std::size_t activeJobs = mAbcCount, pendingJobs = 0;
		while( activeJobs )
		{
			for( auto sidx : waitingSamples )
			{
				auto& sample = mSamples[sidx];

				if( sample.distBis > mGamma && detail::mean_<Scalar>( trials.begin(), trials.end() ) < mAvgTrialCount )
				{
					//XXX-TODO: move this into a separate function (sample_update_()?)
					auto const idx = detail::wrand_idx_( aRng, mPreW.begin(), mPreW.end() );
					assert( idx < mAbcCount );

					for( std::size_t i = 0; i < mComponentCount; ++i )
					{
						sample.mBis[i] = detail::displace_( aRng, mM(i,idx), mTauM[i], mLowerM, mUpperM );
						sample.cBis[i] = detail::displace_( aRng, mC(i,idx), mTauC[i], mLowerC, mUpperC );
					}
					for( std::size_t i = 0; i < mZCount; ++i )
					{
						sample.azBis[i] = detail::displace_( aRng, mAz(i,idx), mTauAz[i], mLowerAz, mUpperAz );
					}

					//XXX- could be templated on ComponentCount
					std::iota( sample.indices.begin(), sample.indices.end(), 0 );
					std::sort( sample.indices.begin(), sample.indices.end(),
						[&sample] (std::size_t aX, std::size_t aY) {
							return sample.mBis[aX] < sample.mBis[aY];
						}
					);

					//XXX- could be templated on ComponentCount
					for( std::size_t i = 0; i < mComponentCount; ++i )
					{
						sample.mTmp[i] = sample.mBis[sample.indices[i]];
						sample.cTmp[i] = sample.cBis[sample.indices[i]];
					}

					std::swap( sample.mTmp, sample.mBis );
					std::swap( sample.cTmp, sample.cBis );

					//XXX- could be templated on ZCount
					if( mZCount > 1 )
					{
						assert( mZCount == mComponentCount );
						for( std::size_t i = 0; i < mComponentCount; ++i )
							sample.azTmp[i] = sample.azBis[sample.indices[i]];

						std::swap( sample.azTmp, sample.azBis );
					}


					Scalar const csum = std::accumulate( sample.cBis.begin(), sample.cBis.end(), Scalar(0) );
					Scalar const lambda = csum * mVolumeFactor;
					std::poisson_distribution<Count> poisson(lambda);

					infoLambdaAcc += lambda;
					infoLambdaCount += Scalar(1);

					assert( sample.particleCounts.size() == mSystemSetup.jobCount );
					for( auto& count : sample.particleCounts )
						count = poisson(aRng);


					std::copy_n( sample.azBis.data(), mZCount, sample.halfAz );

					Scalar acc = Scalar(0);
					for( std::size_t i = 0; i < mComponentCount; ++i )
					{
						acc += sample.cBis[i] / csum;
						sample.preCompProb[i] = acc;
					}

					for( std::size_t i = 0; i < mComponentCount; ++i )
					{
						sample.randWalkStddev[i] = std::sqrt( Scalar(2)*sample.mBis[i]*mDeltaT );
					}

					job_queue_( sidx, sample );
					++pendingJobs;
				}
				else
				{
					assert( activeJobs > 0 );
					--activeJobs;
				}
			}

			waitingSamples.clear();

			if( pendingJobs )
			{
				auto const results = mResults.wait();
				for( auto it = std::get<0>(results); it != std::get<1>(results); ++it )
				{
					assert( it->sample );
					auto& sample = *it->sample;
					auto& dev = mDevGlobal[sample.device];

					CUDA_CHECKED it->error;

					sample.distBis = std::log10( *sample.hostResultDistance );

					sample_dev_clean_( sample, dev );

					auto const sidx = &sample - mSamples.data();
					++trials[sidx];

					--pendingJobs;
					++mInfoSimulations;

					waitingSamples.push_back( sidx );
				}
			}
		}

		infoIterEnd = Clock_::now();

		assert( 0 == activeJobs );
		assert( mResults.empty() );

		if( mVerbosity >= 0 )
		{
			using Fms_ = std::chrono::duration<float,std::milli>;
			
			Count totalTrials = 0;
			Count maxTrials = 0, minTrials = std::numeric_limits<Count>::max();
			
			for( auto x : trials )
			{
				totalTrials += x;
				if( x > maxTrials ) maxTrials = x;
				if( x < minTrials ) minTrials = x;
			}

			std::printf( "Trials: total[min/max]: %4u[%4u/%4u]; γ = %.2g, λ̅ = %.1f\n", totalTrials, minTrials, maxTrials, mGamma, infoLambdaAcc/infoLambdaCount );

			std::printf( "  Wall time: %6.2f ms for this iteration\n", std::chrono::duration_cast<Fms_>(infoIterEnd-infoIterStart).count() );

#			if SIM_KERNEL_TIMINGS
			std::printf( "  GPU: average times: %4.2f ms total (%4.2f sim, %4.2f dist)\n", timeTotalTotal/timeTotalCount, timeSimTotal/timeSimCount, timeDistTotal/timeDistCount );
#			endif // ~ SIM_KERNEL_TIMINGS
		}

		// evaluate weighting scheme; this updates mW
		(this->*mWeightingSchemeFn)();

		// update other state
		for( std::size_t i = 0; i < mAbcCount; ++i )
		{
			auto const& sample = mSamples[i];

			for( std::size_t j = 0; j < mComponentCount; ++j )
			{
				mM(j,i) = sample.mBis[j];
				mC(j,i) = sample.cBis[j];
			}
			for( std::size_t j = 0; j < mZCount; ++j )
			{
				mAz(j,i) = sample.azBis[j];
			}

			mDist[i] = sample.distBis;
		}

		compute_tau_();

		// finished?
		if( detail::mean_<Scalar>( trials.begin(), trials.end() ) >= mAvgTrialCount )
		{
			mConverged = true;
		}
	}
}
template< typename... tArgs > inline
output::Output SimulationT<tArgs...>::output()
{
	output::Output ret;
	ret.model             = kModel;
	ret.weightingScheme   = [&] () {
		if( &SimulationT::weighting_scheme_inv_dist_sq_ == mWeightingSchemeFn )
			return EWeightingScheme::inverseDistSq;
		else if( &SimulationT::weighting_scheme_pmc_standard_ == mWeightingSchemeFn )
			return EWeightingScheme::pmcStandard;

		return EWeightingScheme(~0u);
	}();

	ret.abcSampleCount    = std::size_t(mAbcCount);
	ret.componentCount    = std::size_t(mComponentCount);
	ret.zComponentCount   = std::size_t(mZCount);

	ret.converged         = mConverged;
	ret.gamma             = mGamma;

	ret.m = detail::mat2vec_<double>( mM );
	ret.c = detail::mat2vec_<double>( mC );
	ret.az = detail::mat2vec_<double>( mAz );
	ret.dist.assign( mDist.begin(), mDist.end() );
	ret.w.assign( mW.begin(), mW.end() );

	ret.meta["simulator"] = "gpu";
	ret.meta["iterations"] = tfm::format( "%u", mInfoIterations );
	ret.meta["simulations"] = tfm::format( "%u", mInfoSimulations );
	ret.meta["gpus"] = tfm::format( "%u", mDevGlobal.size() );

	for( std::size_t i = 0; i < mDevGlobal.size(); ++i )
	{
		auto const& dev = mDevGlobal[i];
		ret.meta[tfm::format( "gpu%d_id", i )] = tfm::format( "%d", dev.device );
		ret.meta[tfm::format( "gpu%d_queues", i )] = tfm::format( "%d", dev.queueCount );
	}

	return ret;
}

template< class... tArgs > inline
void SimulationT<tArgs...>::prepare_( input::Parameters const& aPar, HostRng& aRng, SimulationConfig const& aCfg )
{
	if( mComponentCount != aPar.componentCount )
		detail::throw_logic_error_( "SimulationT: component count doesn't match the component count of the input parameters" );
		
	if( kModel != aPar.model )
		detail::throw_logic_error_( "SimulationT: fixed model parameter doesn't match model specified by the input parameters" );
		
	// initialize parameters
	auto const& frameCounts = aPar.frameCounts;
	mKMax = 1 + *std::max_element(frameCounts.begin(), frameCounts.end() );

	mDEBinCount = aPar.deBinCount;
	mDEBinWidth = Scalar(aPar.de.upper / aPar.deBinCount);

	mVolumeFactor = Scalar(aPar.Lx*aPar.Ly*aPar.Lz / Scalar(1e12));

	mLowerM = Scalar(aPar.m.lower);    mUpperM = Scalar(aPar.m.upper);
	mLowerC = Scalar(aPar.c.lower);    mUpperC = Scalar(aPar.c.upper);
	mLowerAz = Scalar(aPar.az.lower);  mUpperAz = Scalar(aPar.az.upper);


	std::uniform_real_distribution<Scalar> mDist( mLowerM, mUpperM );
	std::generate( mM.lbegin(), mM.lend(), [&] { return mDist(aRng); } );

	std::uniform_real_distribution<Scalar> cDist( mLowerC, mUpperC );
	std::generate( mC.lbegin(), mC.lend(), [&] { return cDist(aRng); } );

	std::uniform_real_distribution<Scalar> azDist( mLowerAz, mUpperAz );
	std::generate( mAz.lbegin(), mAz.lend(), [&] { return azDist(aRng); } );


	std::fill( mW.begin(), mW.end(), Scalar(1)/mAbcCount );
	std::partial_sum( mW.begin(), mW.end(), mPreW.begin() );

	compute_tau_();

	mGamma = aPar.adaptiveGamma ? Scalar(-1) : aPar.gamma;
	mDeltaT = Scalar(aPar.deltaT);
	mDeltaGamma = aPar.deltaGamma;
	mAvgTrialCount = Scalar(aPar.avgTrialCount);

	mMaxIter = aCfg.maxIter;
	mVerbosity = aCfg.verbosity;

	// weighting scheme
	mWeightingSchemeFn = nullptr;

	switch( aPar.weightingScheme )
	{
		case EWeightingScheme::pmcStandard: {
			mWeightingSchemeFn = &SimulationT::weighting_scheme_pmc_standard_;
			break;
		}
		case EWeightingScheme::inverseDistSq: {
			mWeightingSchemeFn = &SimulationT::weighting_scheme_inv_dist_sq_;
			break;
		}
	}

	// initialize system settings (passed to kernel)
	mSystemSetup.jobCount = int_cast<Count>(aPar.frameCounts.size());
	mSystemSetup.kmin = int_cast<Count>(aPar.kmin);

	mSystemSetup.deltaT = Scalar(aPar.deltaT);

	mSystemSetup.halfLx = 0.5f * Scalar(aPar.Lx);
	mSystemSetup.halfLy = 0.5f * Scalar(aPar.Ly);
	mSystemSetup.halfLz = 0.5f * Scalar(aPar.Lz);

	mSystemSetup.halfAx = 0.5f * Scalar(aPar.ax);
	mSystemSetup.halfAy = 0.5f * Scalar(aPar.ay);

	// initialize sample
	for( auto& sample : mSamples )
	{
		sample.parent = this;

		resize( sample.mBis, mComponentCount );
		resize( sample.cBis, mComponentCount );
		resize( sample.azBis, mZCount );

		resize( sample.indices, mComponentCount );

		resize( sample.mTmp, mComponentCount );
		resize( sample.cTmp, mComponentCount );
		resize( sample.azTmp, mZCount );

		auto& rd = sample.sampleRunData;
		sample.halfAz          = acquire_host_ptr( rd.halfAz, mZCount );
		sample.preCompProb     = acquire_host_ptr( rd.preCompProb, mComponentCount );
		sample.randWalkStddev  = acquire_host_ptr( rd.randWalkStddev, mComponentCount );

		sample.particleCounts.resize( aPar.frameCounts.size() );
	}

	// generate and upload reference histogram
	Matrix<Count,aspect::MatrixRowMajor> reference( mKMax, mDEBinCount, Count(0) );

	assert( aPar.Ks.size() == aPar.DEs.size() );
	for( std::size_t i = 0; i < aPar.Ks.size(); ++i )
	{
		if( aPar.DEs[i] >= aPar.de.upper )
			continue;

		assert( aPar.DEs[i] >= 0 );
		auto idx = long(aPar.DEs[i] / mDEBinWidth);

		++reference(aPar.Ks[i], idx);
	}

	for( std::size_t c = 0; c < reference.m(); ++c )
	{
		auto cv = reference.col(c);
		std::partial_sum( cv.begin(), cv.end(), cv.begin() );
	}
	for( std::size_t r = 0; r < reference.n(); ++r )
	{
		auto rv = reference.row(r);
		std::partial_sum( rv.begin(), rv.end(), rv.begin() );
	}

	// allocate per-GPU data
	unsigned maxQueueCount = 0;

	char const* spec = aCfg.gpuSpec.c_str();
	while( spec )
	{
		unsigned devID, queueCount = 1;
		int iret = std::sscanf( spec, "%u/%u", &devID, &queueCount );
		if( 1 != iret && 2 != iret )
			detail::throw_invalid_gpuspec_( "Don't understand provided gpuspec" );

		if( mVerbosity >= 2 )
			std::printf( "Note: using device %u with %u queues\n", devID, queueCount );

		--devID;

		if( queueCount > maxQueueCount )
			maxQueueCount = queueCount;

		auto const jobCount = aPar.frameCounts.size();
		std::vector<Count> counts( jobCount );
		for( std::size_t i = 0; i < jobCount; ++i )
			counts[i] = int_cast<Count>(aPar.frameCounts[i]);

		auto const b = (jobCount+4-1)/4;
		auto const t = 32*4;

		Count* fc = nullptr;
		CUDA_CHECKED cudaSetDevice( devID );
		CUDA_CHECKED cudaMalloc( &fc, sizeof(Count)*jobCount );
		CUDA_CHECKED cudaMemcpy( fc, counts.data(), sizeof(Count)*jobCount, cudaMemcpyHostToDevice );

		mDevGlobal.emplace_back( CudaDevGlobal_{
			int(devID),
			queueCount,
			cusim::Histogram2D<Count>( mKMax, mDEBinCount, reference.data() ),
			Pool<Count>( jobCount ),
			MappedPool<Scalar>( 1 ),
			Pool<Scalar>( mZCount ),
			Pool<Scalar>( mComponentCount ),
			Pool<Scalar>( mComponentCount ),
			dim3( b, 1, 1 ),
			dim3( t, 1, 1 ),
			b * t,
			fc
		} );

		if( auto x = std::strchr( spec, ',' ) )
			spec = x+1;
		else
			spec = nullptr;
	}

	if( mVerbosity >= 2 )
		std::printf( "Note: using a total of %zu GPUs\n", mDevGlobal.size() );

	// allocate per-queue data
	for( std::size_t i = 0; i < maxQueueCount; ++i )
	{
		for( std::size_t j = 0; j < mDevGlobal.size(); ++j )
		{
			auto const& devGlobal = mDevGlobal[j];

			if( i >= devGlobal.queueCount )
				continue;

			cudaStream_t stream;
			CUDA_CHECKED cudaSetDevice( devGlobal.device );
			CUDA_CHECKED cudaStreamCreate( &stream );

			mDevQueues.emplace_back( CudaDevQueue_{
				j,
				stream,
				cusim::Histogram2D<Count>( mKMax, mDEBinCount ),
				DeviceRandom::initialize( devGlobal.totalThreads, aRng ),
			} );

			mDevQueues.back().result.clear();
		}
	}

	if( mVerbosity >= 2 )
		std::printf( "Note: using a total of %zu queues\n", mDevQueues.size() );

	if( mDevQueues.empty() )
		detail::throw_invalid_gpuspec_( "No GPU queues were created! Check gpuspec." );
}

template< typename... tArgs > inline
void SimulationT<tArgs...>::compute_tau_()
{
	for( std::size_t i = 0; i < mComponentCount; ++i )
	{
		auto const mr = mM.row(i);
		mTauM[i] = std::sqrt( 2 * detail::var_(mr.begin(), mr.end()) );

		auto const mc = mC.row(i);
		mTauC[i] = std::sqrt( 2 * detail::var_(mc.begin(), mc.end()) );
	}
	for( std::size_t i = 0; i < mZCount; ++i )
	{
		auto const maz = mAz.row(i);
		mTauAz[i] = std::sqrt( 2 * detail::var_(maz.begin(), maz.end()) );
	}
}
template< typename... tArgs > inline
auto SimulationT<tArgs...>::compute_initial_gamma_( HostRng& aRng ) -> Scalar
{
	// XXX-FIXME-quick hack with a supersized side of copy-pasta
	auto copyOfSamples = mSamples;
	std::size_t pendingJobs = 0;

	assert( copyOfSamples.size() == mAbcCount );
	for( std::size_t sidx = 0; sidx < mAbcCount; ++sidx )
	{
		auto& sample = copyOfSamples[sidx];

		for( std::size_t i = 0; i < mComponentCount; ++i )
		{
			sample.mBis[i] = mM(i,sidx);
			sample.cBis[i] = mC(i,sidx);
		}
		for( std::size_t i = 0; i < mZCount; ++i )
		{
			sample.azBis[i] = mAz(i,sidx);
		}

		Scalar const csum = std::accumulate( sample.cBis.begin(), sample.cBis.end(), Scalar(0) );
		Scalar const lambda = csum * mVolumeFactor;
		std::poisson_distribution<Count> poisson(lambda);

		assert( sample.particleCounts.size() == mSystemSetup.jobCount );
		for( auto& count : sample.particleCounts )
			count = poisson(aRng);


		std::copy_n( sample.azBis.data(), mZCount, sample.halfAz );

		Scalar acc = Scalar(0);
		for( std::size_t i = 0; i < mComponentCount; ++i )
		{
			acc += sample.cBis[i] / csum;
			sample.preCompProb[i] = acc;
		}

		for( std::size_t i = 0; i < mComponentCount; ++i )
		{
			sample.randWalkStddev[i] = std::sqrt( Scalar(2)*sample.mBis[i]*mDeltaT );
		}

		job_queue_( sidx, sample );
		++pendingJobs;
	}

	while( pendingJobs )
	{
		auto const results = mResults.wait();
		for( auto it = std::get<0>(results); it != std::get<1>(results); ++it )
		{
			assert( it->sample );
			auto& sample = *it->sample;
			auto& dev = mDevGlobal[sample.device];

			CUDA_CHECKED it->error;

			sample.distBis = std::log10( *sample.hostResultDistance );

			sample_dev_clean_( sample, dev );

			--pendingJobs;
		}
	}

	// find median (log-10) distance
	AVec_<Scalar> distances;
	resize( distances, mAbcCount );

	for( std::size_t i = 0; i < mAbcCount; ++i )
		distances[i] = copyOfSamples[i].distBis;

	std::size_t const n = mAbcCount / 2;
	std::nth_element( distances.begin(), distances.begin()+n, distances.end() );

	return distances[n];
}

template< typename... tArgs > inline
void SimulationT<tArgs...>::weighting_scheme_pmc_standard_()
{
	Scalar wsum = Scalar(0);
	for( std::size_t abc = 0; abc < mAbcCount; ++abc )
	{
		mWStar[abc] = Scalar(0);
		for( std::size_t i = 0; i < mAbcCount; ++i )
		{
			Scalar term = Scalar(1);
			for( std::size_t c = 0; c < mComponentCount; ++c )
			{
				term *= detail::normpdf_( mSamples[abc].mBis[c] - mM(c,i), Scalar(0), mTauM[c] );
				term *= detail::normpdf_( mSamples[abc].cBis[c] - mC(c,i), Scalar(0), mTauC[c] );
			}
			for( std::size_t z = 0; z < mZCount; ++z )
			{
				term *= detail::normpdf_( mSamples[abc].azBis[z] - mAz(z,i), Scalar(0), mTauAz[z] );
			}

			mWStar[abc] += mW[abc] * term;
		}

		mWStar[abc] = Scalar(1) / mWStar[abc];
		wsum += mWStar[abc];
	}

	for( std::size_t abc = 0; abc < mAbcCount; ++abc )
		mW[abc] = mWStar[abc] / wsum;
}
template< typename... tArgs > inline
void SimulationT<tArgs...>::weighting_scheme_inv_dist_sq_()
{
	Scalar wsum = Scalar(0);
	for( std::size_t i = 0; i < mAbcCount; ++i )
	{
		auto const& sample = mSamples[i];
		auto const w = Scalar(1) / (sample.distBis*sample.distBis);

		mW[i] = w;
		wsum += w;
	}

	for( auto& w : mW )
		w /= wsum;
}

template< typename... tArgs > inline
void SimulationT<tArgs...>::job_queue_( std::size_t, Sample_& aSample )
{
	//XXX-round robin. perhaps do something more dynamic...
	auto const queueIdx = mNextQueue++;
	if( mNextQueue >= mDevQueues.size() )
		mNextQueue = 0;

	auto& queue = mDevQueues[queueIdx];
	auto& dev = mDevGlobal[queue.devidx];

	aSample.device = queue.devidx;
	sample_dev_init_( aSample, dev, queue );

	CUDA_CHECKED cudaSetDevice( dev.device );

	{
#		if SIM_KERNEL_TIMINGS
		cudaEventRecord( aSample.simStart, queue.stream );
#		endif // ~ SIM_KERNEL_TIMINGS

		cusim::K_simulate_system<<<dev.blocks,dev.threads,0,queue.stream>>>(
			mSystemSetup,
			aSample.sampleRunData,
			cusim::HistogramRecorder<Count,Scalar>(queue.result,mDEBinWidth),
			queue.randomState
		);

#		if SIM_KERNEL_TIMINGS
		cudaEventRecord( aSample.simStop, queue.stream );
#		endif // ~ SIM_KERNEL_TIMINGS
	}

	{
#		if SIM_KERNEL_TIMINGS
		cudaEventRecord( aSample.distStart, queue.stream );
#		endif // ~ SIM_KERNEL_TIMINGS

		dim3 const bl(1,1,1), th(32,32,1);
		cusim::K_distance<<<bl,th,0,queue.stream>>>(
			aSample.devResultDistance,
			int_cast<unsigned>(mKMax),
			int_cast<unsigned>(mDEBinCount),
			queue.result.device_ptr(),
			dev.reference.device_ptr()
		);

#		if SIM_KERNEL_TIMINGS
		cudaEventRecord( aSample.distStop, queue.stream );
#		endif // ~ SIM_KERNEL_TIMINGS
	}

#	if SIM_KERNEL_TIMINGS
	cudaEventRecord( aSample.totalStop, queue.stream ); // WARN: stared by dev_init_()
#	endif // ~ SIM_KERNEL_TIMINGS

	CUDA_CHECKED cudaStreamAddCallback(
		queue.stream,
		&SimulationT::cuda_stream_callback_,
		const_cast<Sample_*>(&aSample),
		0
	);
}

template< typename... tArgs > inline
void SimulationT<tArgs...>::sample_dev_init_( Sample_& aSample, CudaDevGlobal_& aDev, CudaDevQueue_& aQueue )
{
	CUDA_CHECKED cudaSetDevice( aDev.device );

#	if SIM_KERNEL_TIMINGS
	aSample.simStart = aDev.timeEvents.alloc();
	aSample.simStop = aDev.timeEvents.alloc();
	aSample.distStart = aDev.timeEvents.alloc();
	aSample.distStop = aDev.timeEvents.alloc();
	aSample.totalStart = aDev.timeEvents.alloc();
	aSample.totalStop = aDev.timeEvents.alloc();

	cudaEventRecord( aSample.totalStart, aQueue.stream ); // WARN: must end outside.
#	endif // ~ SIM_KERNEL_TIMINGS

	// NOTE: the histogram is cleared on allocation, and K_distance zeroes it
	// before returning -- there's no need to clear it one more time.
	//aQueue.result.clear_async( aQueue.stream );

	upload_from_host_ptr(
		aSample.sampleRunData.halfAz,
		aSample.halfAz,
		mZCount,
		[&aDev] () { return aDev.halfAzPool.alloc(); },
		aQueue.stream
	);
	upload_from_host_ptr(
		aSample.sampleRunData.preCompProb,
		aSample.preCompProb,
		mComponentCount,
		[&aDev] () { return aDev.preCompProbPool.alloc(); },
		aQueue.stream
	);
	upload_from_host_ptr(
		aSample.sampleRunData.randWalkStddev,
		aSample.randWalkStddev,
		mComponentCount,
		[&aDev] () { return aDev.randWalkStddevPool.alloc(); },
		aQueue.stream
	);

	auto const particleCounts = aDev.particleCountPool.alloc();
	CUDA_CHECKED cudaMemcpyAsync(
		particleCounts,
		aSample.particleCounts.data(),
		sizeof(Count)*aSample.particleCounts.size(),
		cudaMemcpyHostToDevice,
		aQueue.stream
	);

	auto hp = aDev.resultDistancePool.alloc();
	aSample.hostResultDistance = std::get<0>(hp);
	aSample.devResultDistance = std::get<1>(hp);

	aSample.sampleRunData.frames = aDev.devFrameCounts;
	aSample.sampleRunData.particles = particleCounts;
}
template< typename... tArgs > inline
void SimulationT<tArgs...>::sample_dev_clean_( Sample_& aSample, CudaDevGlobal_& aDev )
{
	CUDA_CHECKED cudaSetDevice( aDev.device );

	aDev.particleCountPool.free( aSample.sampleRunData.particles );
	aDev.resultDistancePool.free( aSample.hostResultDistance, aSample.devResultDistance );

	clean_gpu_cache(
		aSample.sampleRunData.halfAz,
		[&aDev] (Scalar* aPtr) { aDev.halfAzPool.free( aPtr ); }
	);
	clean_gpu_cache(
		aSample.sampleRunData.preCompProb,
		[&aDev] (Scalar* aPtr) { aDev.preCompProbPool.free( aPtr ); }
	);
	clean_gpu_cache(
		aSample.sampleRunData.randWalkStddev,
		[&aDev] (Scalar* aPtr) { aDev.randWalkStddevPool.free( aPtr ); }
	);

#	if SIM_KERNEL_TIMINGS
	{
		float sim, dist, tot;
		CUDA_CHECKED cudaEventElapsedTime( &sim, aSample.simStart, aSample.simStop );
		CUDA_CHECKED cudaEventElapsedTime( &dist, aSample.distStart, aSample.distStop );
		CUDA_CHECKED cudaEventElapsedTime( &tot, aSample.totalStart, aSample.totalStop );

		timeSimTotal += sim; timeSimCount += 1;
		timeDistTotal += dist; timeDistCount += 1;
		timeTotalTotal += tot; timeTotalCount += 1;

		aDev.timeEvents.free( aSample.simStart );
		aDev.timeEvents.free( aSample.simStop );
		aDev.timeEvents.free( aSample.distStart );
		aDev.timeEvents.free( aSample.distStop );
		aDev.timeEvents.free( aSample.totalStart );
		aDev.timeEvents.free( aSample.totalStop );
	}
#	endif // ~ SIM_KERNEL_TIMINGS
}

template< typename... tArgs > inline
void SimulationT<tArgs...>::cuda_stream_callback_( cudaStream_t, cudaError_t aError, void* aUser )
{
	/* NOTE: callback from CUDA. We're *not* allowed to perform any CUDA calls
	 * in here! Furthermore, this may be called from a thread different from the
	 * main program.
	 */

	auto& sample = *static_cast<Sample_*>(aUser);
	auto* self = sample.parent;

	self->mResults.queue( Result_{
		&sample,
		aError
	} );
}
