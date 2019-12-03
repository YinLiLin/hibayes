#ifndef MYTIMER_H
#define MYTIMER_H
#include <stdint.h>
#include <vector>
#include <string>
#define R_NO_REMAP
#include <Rinternals.h>
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#elif defined(__APPLE__)
#include <mach/mach_time.h>
#elif defined(linux) || defined(__linux) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__GLIBC__) || defined(__GNU__) || defined(__CYGWIN__)
#include <time.h>
#elif defined(sun) || defined(__sun) || defined(_AIX)
#include <sys/time.h>
#else /* Unsupported OS */
#error "Rcpp::MyTimer not supported by your OS."
#endif
namespace Rcpp{
    typedef uint64_t nanotime_t;
	#if defined(_WIN32)
    inline nanotime_t get_nanotime(void) {
        LARGE_INTEGER time_var, frequency;
        QueryPerformanceCounter(&time_var);
        QueryPerformanceFrequency(&frequency);
        return 1.0e9 * time_var.QuadPart / frequency.QuadPart;
    }
	#elif defined(__APPLE__)
    inline nanotime_t get_nanotime(void) {
        nanotime_t time;
        mach_timebase_info_data_t info;
        time = mach_absolute_time();
        mach_timebase_info(&info);
        return time * (info.numer / info.denom);
    }
	#elif defined(linux) || defined(__linux) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__GLIBC__) || defined(__GNU__) || defined(__CYGWIN__)
    static const nanotime_t nanoseconds_in_second = static_cast<nanotime_t>(1000000000.0);
    inline nanotime_t get_nanotime(void) {
        struct timespec time_var;
        clock_gettime(CLOCK_REALTIME, &time_var);
        nanotime_t sec = time_var.tv_sec;
        nanotime_t nsec = time_var.tv_nsec;
		return (nanoseconds_in_second * sec) + nsec;
    }
	#elif defined(sun) || defined(__sun) || defined(_AIX)
    inline nanotime_t get_nanotime(void) {
        return gethrtime();
    }

	#endif
    class MyTimer {
    public:
        MyTimer() : data(2), start_time( get_nanotime() ){}
        MyTimer(nanotime_t start_time_) : data(), start_time(start_time_){}
        void step( const std::string& name){
            if(name == "start"){data[0] = (std::make_pair(name, now()));}
            if(name == "end"){data[1] = (std::make_pair(name, now()));}
        }
        operator SEXP() const {
            size_t n = data.size();
            NumericVector out(n);
            CharacterVector names(n);
            for (size_t i=0; i<n; i++) {
                names[i] = data[i].first;
                out[i] = data[i].second - start_time ;
            }
            out.attr("names") = names;
            return out;
        }
        static std::vector<MyTimer> get_timers(int n){
            return std::vector<MyTimer>( n, MyTimer() ) ;
        }
        inline nanotime_t now() const {
            return get_nanotime() ;
        }
        inline nanotime_t origin() const {
            return start_time ;
        }
    	private:
        typedef std::pair<std::string,nanotime_t> Step;
        typedef std::vector<Step> Steps;
        Steps data;
        const nanotime_t start_time;
    };
}
#ifdef FALSE
#undef FALSE
#endif
#endif
