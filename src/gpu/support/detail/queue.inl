#include <utility>

template< typename tType > template< typename... tArgs > inline
void QueueSink<tType>::queue( tArgs&&... aArgs )
{
	{
		std::unique_lock<std::mutex> lock(mMutex);
		mLive.emplace_back( std::forward<tArgs>(aArgs)... );
	}
	
	mCV.notify_one();
}

template< typename tType > inline
bool QueueSink<tType>::empty() const
{
	std::unique_lock<std::mutex> lock(mMutex);
	return mLive.empty();
}



template< typename tType > inline
std::tuple<tType*,tType*> Queue<tType>::wait()
{
	this->mDead.clear();
	
	{
		std::unique_lock<std::mutex> lock(this->mMutex);
		this->mCV.wait( lock, [this] { return !this->mLive.empty(); } );

		std::swap( this->mLive, this->mDead );
	}

	auto ptr = this->mDead.data();
	return std::make_tuple( ptr, ptr+this->mDead.size() );
}
template< typename tType > inline
std::tuple<tType*,tType*> Queue<tType>::poll()
{
	this->mDead.clear();
	
	{
		std::unique_lock<std::mutex> lock(this->mMutex);
		std::swap( this->mLive, this->mDead );
	}

	auto ptr = this->mDead.data();
	return std::make_tuple( ptr, ptr+this->mDead.size() );
}

