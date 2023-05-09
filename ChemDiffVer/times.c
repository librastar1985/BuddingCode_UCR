/*		
*			times.c
*
*			Routines for timing purposes.
*
*	For the most part we define the cpu time of a process as the
*	sum of its user and system times and of those of all its terminated
*	subprocesses.
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include "cdecs.h"

#include <sys/types.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

#if !defined(CLK_TCK)
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif /* !defined(CLK_TCK) */

/* LINTLIBRARY */

	/* LOCAL Function prototypes*/
LOCAL	void	print_rss_usage(void);

LOCAL	int  MAX_TIMES = 0;	/* Length of cputime stack */
LOCAL	float *cputime = NULL;	/* Stack for  Storing Recursive Times */
LOCAL	float *walltime = NULL;	/* Stack for  Storing Recursive Times */
LOCAL int  top = 0;			/* Pointer To top of Stack */

EXPORT void print_execution_times(void)
{
	struct tms Tm;
	size_t	   usert_min, usert_sec, usert_day, usert_hr;
	size_t	   syst_min, syst_sec, syst_day, syst_hr;
	size_t     clk_tck;

	clk_tck = CLK_TCK;
	(void) times(&Tm);

	usert_sec = (Tm.tms_utime+Tm.tms_cutime)/clk_tck;
	usert_day = (usert_sec/86400);
	usert_hr = (usert_sec/3600) % 24;
	usert_min = (usert_sec/60) % 60;
	usert_sec = usert_sec%60;
	syst_sec = (Tm.tms_stime+Tm.tms_cstime)/clk_tck;
	syst_day = (syst_sec/86400);
	syst_hr = (syst_sec/3600) % 24;
	syst_min = (syst_sec/60) % 60;
	syst_sec = syst_sec%60;

	screen("\n\n\t\tEXECUTION TIME DATA:\n\n\t\tCPU TIME: ");
	if (usert_day) screen("%d Day%s ",usert_day,(usert_day==1)?"":"s");
	if (usert_hr) screen("%d Hour%s ",usert_hr,(usert_hr==1)?"":"s");
	if (usert_min) screen("%d Minute%s ",usert_min,(usert_min==1)?"":"s");
	screen("%d Second%s",usert_sec,(usert_sec==1)?"":"s");
	screen("\n\t\tSYS TIME: ");
	if (syst_day) screen("%d Day%s ",syst_day,(syst_day==1)?"":"s");
	if (syst_hr) screen("%d Hour%s ",syst_hr,(syst_hr==1)?"":"s");
	if (syst_min) screen("%d Minute%s ",syst_min,(syst_min==1)?"":"s");
	screen("%d Second%s\n",syst_sec,(syst_sec==1)?"":"s");
	print_rss_usage();
}		/*end print_execution_times*/





EXPORT void start_clock(
	const char	*s)
{
	int	i;

	if (!debugging("CLOCK"))
	    return;

	if (cputime == NULL)
	{
	    MAX_TIMES = 40;
	    cputime = (float*)malloc(MAX_TIMES*sizeof(float));
	    walltime = (float*)malloc(MAX_TIMES*sizeof(float));	    
	    zero_scalar(cputime,MAX_TIMES*sizeof(float));
	    zero_scalar(walltime,MAX_TIMES*sizeof(float));
	}

	if  (top >= MAX_TIMES)
	{
	    float	*new_cputime;
	    float	*new_walltime;	    
	    int	NEW_MAX_TIMES = 2*MAX_TIMES;

	    new_cputime = (float*)malloc(NEW_MAX_TIMES*sizeof(float));
	    new_walltime = (float*)malloc(NEW_MAX_TIMES*sizeof(float));	    
	    zero_scalar(new_cputime,NEW_MAX_TIMES*sizeof(float));
	    zero_scalar(new_walltime,NEW_MAX_TIMES*sizeof(float));	    
	    for (i = 0; i < MAX_TIMES; i++) 
	    	new_cputime[i] = cputime[i];
	    for (i = 0; i < MAX_TIMES; i++) 
	    	new_walltime[i] = walltime[i];	    
	    MAX_TIMES = NEW_MAX_TIMES;
            (void)free(cputime); 
            (void)free(walltime); 
            cputime = new_cputime;
            walltime = new_walltime;  
	}

	for (i = 0; i < top; i++)
	    (void) printf("  ");
	(void) printf("CLOCK           (%s)\n",s);
	cputime[top] = cpu_seconds();
	walltime[top] = real_time();	
	top++;
}		/*end start_clock*/


EXPORT void stop_clock(
	const char	*s)
{
	if (!debugging("CLOCK"))
	    return;
	if (cputime == NULL)
	    return;
	top--;
	if (top < 0)
	{
	    (void) printf("ERROR: stop_clock(%s): CLOCK STACK EMPTY\n",s);
	    top = 0;
	}
	else if (top >= MAX_TIMES) 
	    (void) printf("ERROR: stop_clock(%s): CLOCK STACK FULL\n",s);
	else
	{
	    int i;
	    for (i = 0; i < top; i++)
		(void) printf("  ");
	    (void) printf("CLOCK  %7.3f   (%s) %7.3f\n", /* %6.2f */
	    		  (cpu_seconds() - cputime[top]),s,
	    		  (real_time() - walltime[top]));
	}
}		/*end stop_clock*/




EXPORT void cpu_time(
	const char	*s)
{
	(void) printf("TIME: %7.2f  at %s\n", cpu_seconds(), s);
}		/*end cpu_time*/




/*
*				cpu_seconds():
*
*	Returns the number of CPU seconds used by the program to this point
*	Includes both the user and system time for both the program and
*	all its terminated subprocesses.
*/


EXPORT double cpu_seconds(void)
{
	struct tms	Tm;
	double		cpu_time;
	double		clk_tck;

	(void) times(&Tm);
	clk_tck = (double)CLK_TCK;
	cpu_time = Tm.tms_utime + Tm.tms_cutime;
	return cpu_time/clk_tck;
}		/*end cpu_seconds*/




/*
*				real_time():
*
*	Returns the current real time in seconds measured from
*	some constant date.
*/

EXPORT float real_time(void)
{
        /*
	time_t		tvec;
	static time_t	stvec;
	static int	first = 1;

	if (first)
	{
		first = 0;
		(void) time(&stvec);
	}

	(void) time(&tvec);
	return (float)(tvec-stvec);
        */
        /* The new code gives the fractional seconds
         * up to microsecond */

        static int      first = 1;
        struct timeval  Tval;
        static double   tsec;

        if (first)
        {
            first = 0;
            gettimeofday (&Tval, NULL);
            tsec = (double)Tval.tv_sec + (double)(Tval.tv_usec)/1000000.0;
        }
        gettimeofday (&Tval, NULL);

        return (float)((double)Tval.tv_sec + 
                       (double)(Tval.tv_usec)/1000000.0 - tsec);
}		/*end real_time*/



EXPORT char *date_string(void)
{
	time_t	tvec;

	time (&tvec);
	return ctime(&tvec);
}		/*end date_string*/


#if !(defined(__SUNPRO_C) || defined(__SUNPRO_CC)) && !defined(_HPUX_SOURCE) && !defined(cray) && !defined(tflops)
#include <sys/time.h>
#include <sys/resource.h>

LOCAL	void	print_rss_usage(void)
{
	struct rusage	Rusg, Crusg;
	int		pagesz = 1;

	(void) getrusage(RUSAGE_SELF,&Rusg);
	(void) getrusage(RUSAGE_CHILDREN,&Crusg);

#if defined(sparc)
	pagesz = getpagesize()/1024;
#endif /* defined(sparc) */
	screen("\t\tMAX RSS = %d Kilobytes\n",
		(Rusg.ru_maxrss+Crusg.ru_maxrss)*pagesz);
}		/*end print_rss_usage*/
#else /* !(defined(__SUNPRO_C) || defined(__SUNPRO_CC)) && !defined(_HPUX_SOURCE) && !defined(cray) && !defined(tflops) */

/*ARGSUSED*/
LOCAL	void	print_rss_usage(void)
{
}		/*end print_rss_usage*/

#endif /* !(defined(__SUNPRO_C) || defined(__SUNPRO_CC)) && !defined(_HPUX_SOURCE) && !defined(cray) && !defined(tflops) */

