#ifndef CT_TIMER_H
#define CT_TIMER_H
#include <ctime>
#include <chrono>
#include <string>
#include <boost/format.hpp>

using namespace std::chrono;

namespace Cantera {
//! The class provides a timer in seconds
/*!
 * This routine relies on the chrono library for its basic operation.
 * Therefore, it should be fairly portable.
 *
 *
 * An example of how to use the timer is given below. timeToDoCalcs contains the
 * wall clock time calculated for the operation.
 *
 * @code
 * Timer* timer_test = new Timer("timer_test");
 * timer_test->start();
 * computations()
 * timer_test->stop();
 * timer_test->print();
 * @endcode
 *
 *
 * @ingroup globalUtilFuncs
 *
 */
class Timer
{
public:
    Timer() {}
    Timer(const std::string name) {
        m_name = name;
    }

	//! Start of code block to monitor timer
    void start() {
        m_tstart = steady_clock::now();
        m_count++;
    }

	//! End of code block to monitor timer
    void stop() {
        m_tend = steady_clock::now();
        m_time_span += duration_cast<duration<double>>(m_tend - m_tstart);
    }

    //! Return string containing stats of the Timer object
    std::string stats() {
        std::string stats;
        stats = boost::str(boost::format("---> %s - count = %d - time spent = %.3e s") 
                                % m_name % m_count % m_time_span.count());
        return stats;
    }

    int count() {
        return m_count;
    }
private:
    steady_clock::time_point m_tstart, m_tend;
    duration<double> m_time_span = steady_clock::duration::zero();
    int m_count = 0;
    std::string m_name = "timer";
};
}
#endif