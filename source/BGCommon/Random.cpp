#include <cinttypes>
#include <ctime>
#include <thread>
#include <unistd.h>

#include "Random.hpp"

uint32_t generateSeed()
{
	union
	{
		uint64_t t64;
		uint32_t t32[2];
	};
	t64 = 179424743ll*time(0);
	uint32_t seed = t32[0]^t32[1];
	timespec ts;
	clock_gettime(CLOCK_MONOTONIC,&ts);
	seed ^= 179425153ll*ts.tv_nsec;
	seed ^= 179426341ll*getpid();
	std::thread::id this_id = std::this_thread::get_id();
	seed ^= static_cast<uint32_t>(179426339ll*std::hash<std::thread::id>{}(this_id));
	return seed;
}
