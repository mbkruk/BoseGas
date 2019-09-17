#ifndef BARRIER_HPP_
#define BARRIER_HPP_

#include <cinttypes>
#include <thread>
#include <mutex>
#include <condition_variable>

class Barrier
{
private:

	uint_fast32_t threadCount;
	uint_fast32_t counter;
	uint_fast32_t generation = 0;
	std::mutex m;
	std::condition_variable cv;

public:

	inline void wait()
	{
		std::unique_lock<std::mutex> lk(m);
		uint_fast32_t gen = generation;
		if(!--counter)
		{
			++generation;
			counter = threadCount;
			cv.notify_all();
		}
		else
		{
			cv.wait(lk,[this,gen]{return gen!=generation;});
		}
	}

	inline void operator() ()
	{
		wait();
	}

	inline void setThreadCount(uint_fast32_t cnt)
	{
		threadCount = cnt;
		counter = threadCount;
	}

	inline uint_fast32_t getThreadCount() const
	{
		return threadCount;
	}

	inline uint_fast32_t getGeneration() const
	{
		return generation;
	}

	Barrier(const Barrier&) = delete;

	inline Barrier(uint_fast32_t threadCount_=1)
	{
		setThreadCount(threadCount_);
	}
};

#endif
