#ifndef CLOCK_H
#define CLOCK_H

//clock() function; CLOCKS_PER_SEC constant
#include <time.h>
#include "BasicTypes.h"             //uint, ulong
#include <iostream> //for cout
#ifdef __linux__
#include <sys/time.h> //for timing
#elif _WIN32
#include <time.h>
#endif

struct Clock
{
  Clock(): stepsPerOut(0), prevStep(0), lastStep(0), lastTime(0.0) {}
  void Init(const ulong steps, const ulong totSt)
  {
    stepsPerOut = steps;
#ifdef __linux__
    gettimeofday(&tv, &tz);
    strt = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
#elif _WIN32
    strt = clock();
#endif
    lastTime = strt;
    lastStep = totSt - 1;
  }
  void CheckTime(const uint step);
private:
  double TimeInSec(const double strt, const double stp)
  {
    return (double(stp)-double(strt))/CLOCKS_PER_SEC;
  }

#ifdef __linux__
  struct timeval tv;
  struct timezone tz;
  double strt, stop, lastTime;
#elif _WIN32
  clock_t strt, lastTime;
#endif
  ulong stepsPerOut, prevStep, lastStep;
};

inline void Clock::CheckTime(const uint step)
{
  uint stepDelta = step - prevStep;
  if (stepDelta == stepsPerOut && step != lastStep)
  {
#ifdef __linux__
    gettimeofday(&tv, &tz);
    double currTime = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
    std::cout << "Steps/sec. : "
              << stepDelta/(currTime - lastTime) << std::endl;
#elif _WIN32
    clock_t currTime = clock();
    std::cout << "Steps/sec. : "
              << stepDelta/((double)currTime - lastTime/CLOCKS_PER_SEC)
              << std::endl;
#endif
    prevStep = step;
    lastTime = currTime;
  }
  else if (step == lastStep)
  {
#ifdef __linux__
    gettimeofday(&tv, &tz);
    stop = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
    std::cout << "Simulation Time (total): " << (stop - strt)
              << " sec." << std::endl;
#elif _WIN32
    clock_t stop = clock();
    std::cout << "Simulation Time (total): "
              << ((double)stop - strt/CLOCKS_PER_SEC)
              << " sec." << std::endl;
#endif

  }
}

#endif /*CLOCK_H*/
