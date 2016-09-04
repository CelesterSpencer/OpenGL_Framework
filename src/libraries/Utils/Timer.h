//
// Created by ubundrian on 22.08.16.
//

#ifndef OPENGL_FRAMEWORK_TIMER_H
#define OPENGL_FRAMEWORK_TIMER_H

#include <chrono>
#include <deque>
#include <string>

#include "Logger.h"

class Timer {
public:
    std::chrono::steady_clock::time_point startTime;
    std::chrono::steady_clock::time_point endTime;
    std::deque<long> durations;
    long totalDuration = 0;
    int limit = 100;

    void start();
    void stop();
    double getDuration();
};


#endif //OPENGL_FRAMEWORK_TIMER_H
