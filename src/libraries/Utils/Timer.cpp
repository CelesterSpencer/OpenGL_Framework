//
// Created by ubundrian on 22.08.16.
//

#include "Timer.h"

void Timer::start()
{
    startTime = std::chrono::steady_clock::now();
    endTime = startTime;
}

void Timer::stop()
{
    endTime = std::chrono::steady_clock::now();
    long duration = std::chrono::duration_cast<std::chrono::microseconds>( endTime - startTime ).count();
    durations.push_back(duration);
    totalDuration += duration;
    if (durations.size() >= limit) {
        long oldestDuration = durations.front();
        totalDuration -= oldestDuration;
        durations.pop_front();
    }
}

double Timer::getDuration()
{
    double res = 0;
    for (int i = 0; i < durations.size(); i++) {
        res += (double)durations.at(i);
    }
    return (res/durations.size()) / 1000.0;
    //return (totalDuration / durations.size()) / 1000.0;
}