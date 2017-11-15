/*-******************************************************* -- HEADER -{{{1- */
/*-	Simple multi-threaded producer/consumer queue
 *
 * Targets a use case with multiple producers and a single consumer, where 
 * synchronization is unlikely to be a large overhead (i.e., no attempts at
 * making this lock-free).
 */
/*-***************************************************************** -}}}1- */

#ifndef QUEUE_HPP_9A638130_C91E_4DBC_A9AD_17763309CCDE
#define QUEUE_HPP_9A638130_C91E_4DBC_A9AD_17763309CCDE

#include <tuple>
#include <mutex>
#include <vector>
#include <condition_variable>

template< typename tType >
class QueueSink
{
	public: 
		// Note: it is safe to call queue() from multiple threads and/or
		// concurrently to wait()/poll().
		template< typename... tArgs >
		void queue( tArgs&&... );

		bool empty() const;

	protected:
		std::vector<tType> mLive;
		std::vector<tType> mDead;
	
		mutable std::mutex mMutex;
		std::condition_variable mCV;
};

template< typename tType >
class Queue final : public QueueSink<tType> 
{
	public:
		/* Note: wait() and poll() both invalidate their previous return
		 * values. One must handle the full returned range before calling
		 * either function again.
		 *
		 * WARNING: wait() and poll() have race conditions if called from
		 * multiple threads!
		 */
	
		std::tuple<tType*,tType*> wait();
		std::tuple<tType*,tType*> poll();
};

#include "detail/queue.inl"
#endif // QUEUE_HPP_9A638130_C91E_4DBC_A9AD_17763309CCDE
