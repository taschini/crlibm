#ifndef TBX_TIMING_H
#define TBX_TIMING_H

#include <sys/types.h>
#include <sys/time.h>

/*
 * For the user interested in High resolution timer, to my knowledge 
 * the best library is "RTAI", unfortunately it works only on Linux OS
 */


typedef union u_tbx_tick
{
  unsigned long long tick;

  struct
  {
    unsigned long low;
    unsigned long high;
  } sub;

  struct timeval timev;
} tbx_tick_t, *p_tbx_tick_t;

#if defined(SCS_TYPECPU_ITANIUM)
#define TBX_GET_TICK(t) \
   __asm__ __volatile__("mov %0=ar.itc" : "=r"((t).tick) :: "memory")
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2).tick - (t1).tick)

/* Hum hum... Here we suppose that X86ARCH => Pentium! */
#elif defined(SCS_TYPECPU_X86)
#define TBX_GET_TICK(t) \
   __asm__ volatile("cpuid; rdtsc" : "=a" ((t).sub.low), "=d" ((t).sub.high))
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2).tick - (t1).tick)

#elif defined(SCS_TYPECPU_ALPHA)
#define TBX_GET_TICK(t) \
   __asm__ volatile("rpcc %0\n\t" : "=r"((t).tick))
#define TBX_TICK_RAW_DIFF(t1, t2) \
   (((t2).tick & 0xFFFFFFFF) - ((t1).tick & 0xFFFFFFFF))

#elif defined(SCS_TYPECPU_SPARC)
#define TBX_GET_TICK(t) \
    (t).tick = gethrtime()
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2).tick  - (t1).tick)
/*
The following instructions are done only for ppc 601 (which is an old processor!)
#elif defined(SCS_TYPECPU_POWERPC)
#define TBX_GET_TICK(t) \
   __asm__ volatile("li %0,64;mtspr UMMCR0,%0;mfspr %0,UPMC1": "=r"((t).sub.low))
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2).tick  - (t1).tick)


The following one is extracted from Motorola reference manual for 32 bits PPCs.
But apparently it give randomly numbers, may be the clock generator is a random generator !
Silly processor !!!!

#elif defined(SCS_TYPECPU_POWERPC)
#define TBX_GET_TICK(t) \
 {unsigned long chk; \
   __asm__ volatile("0: mftbu %0; mftb %1; mftbu %2; cmpw %2, %0; bne 0b" \
		    : "=r" ((t).sub.low), "=r" ((t).sub.high), "=r" (chk) );}
*/ 
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2).tick  - (t1).tick)

#else
#define TBX_GET_TICK(t) \
   gettimeofday(&(t).timev, 0)
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2.timev.tv_sec * 1000000L + t2.timev.tv_usec) - \
    (t1.timev.tv_sec * 1000000L + t1.timev.tv_usec))

#endif

#define TBX_TICK_DIFF(t1, t2) (TBX_TICK_RAW_DIFF(t1, t2) - tbx_residual + 1)
#define TBX_TIMING_DELAY(t1, t2) tbx_tick2usec(TBX_TICK_DIFF(t1, t2))

extern unsigned long long tbx_residual;
extern tbx_tick_t         tbx_new_event;
extern tbx_tick_t         tbx_last_event;

#endif /* TBX_TIMING_H */
