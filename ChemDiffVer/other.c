/*
*			other.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include "cdecs.h"
#include <stdarg.h>

enum { MAX_FLOATS = 255 };

/*
*			print_line_of_floats():
*
*	Prints a collection of floats in binary or non-binary format
*	See the comments for print_vector_of_floats().
*
*	usage:  print_line_of_floats(n,f1,f2,f3,...,fn)
*
*	Useful for printing columnar data
*/



/* VARARGS */

EXPORT void print_line_of_floats(
	int	number,
	...)
{
	va_list ap;
	int	num,i;
	double	next;
	float	f[MAX_FLOATS];
	       /* MAX_FLOATS = maximum number of floats in line */

	va_start(ap, number);


	if (is_binary_output() == YES)
	{
	    while (number)
	    {
	    	num = (number>MAX_FLOATS) ? MAX_FLOATS : number;
	    	(void) printf("\f%c",num);
	    	for (i = 0; i < num; ++i)
	    	{
	    	    next = va_arg(ap,double);
	    	    f[i] = (float)next;
	    	}
	    	(void) fwrite((const void *)f,sizeof(float),num,stdout);
	    	number -= num;
	    }
	}
	else
	{
	    for (i = 0; i < number; ++i)
	    {
	    	next = va_arg(ap,double);
	    	(void) printf("%-"FFMT" ",next);
	    }
	    (void) printf("\n");
	}
	va_end(ap);
}		/*end print_line_of_floats*/

/*
*			fprint_line_of_floats():
*
*	Prints a collection of floats in binary or non-binary format
*	See the comments for print_vector_of_floats().
*
*	usage:  fprint_line_of_floats(file,n,f1,f2,f3,...,fn)
*
*	Useful for printing columnar data
*/



/* VARARGS */

EXPORT void fprint_line_of_floats(
	FILE	*file,
	int	number,
	 ...)
{
	va_list		ap;
	int		num,i;
	double		next;
	float		f[MAX_FLOATS];
		       /* MAX_FLOATS = maximum number of floats in line */

	va_start(ap, number);


	if (is_binary_output() == YES)
	{
	    while (number)
	    {
	    	num = (number>MAX_FLOATS) ? MAX_FLOATS : number;
		(void) fprintf(file,"\f%c",num);
		for (i = 0; i < num; ++i)
		{
		    next = va_arg(ap,double);
		    f[i] = (float)next;
		}
		(void) fwrite((const void *)f,sizeof(float),num,file);
		number -= num;
	    }
	}
	else
	{
	    for (i = 0; i < number; ++i)
	    {
	    	next = va_arg(ap,double);
	    	(void) fprintf(file,"%-"FFMT" ",next);
	    }
	    (void) fprintf(file,"\n");
	}
	va_end(ap);
}		/*end fprint_line_of_floats*/

EXPORT	const char *right_flush(
	int	   n,
	int	   ndigits)
{
        static  char    s[20], tmp[20];
        int             i;

        if (n == 0)
                ndigits--;
        for (i = n; i > 0; i /= 10) ndigits--;

        s[0] = '\0';
        for (;ndigits > 0; ndigits--)
        {
            (void) sprintf(tmp,"%s0",s);
            (void) strcpy(s,tmp);
        }
        (void) sprintf(tmp,"%s%d",s,n);
        (void) strcpy(s,tmp);
        /*
        for (;ndigits > 0; ndigits--)
                (void) sprintf(s,"%s0",s);
        (void) sprintf(s,"%s%d",s,n);
        */
        return s;
}		/*end right_flush*/

EXPORT	const	char	*y_or_n(
	bool	arg)
{
	switch (arg)
	{
	case YES:
	    return "yes";
	case NO:
	default:
	    return "no";
	}
}		/*end y_or_n*/

EXPORT	const char *ordinal_suffix(
	int i)
{
	i = (i >= 0) ? i%20 : (-i)%20;
	switch (i)
	{
	case 1:
	    return "st";
	case 2:
	    return "nd";
	case 3:
	    return "rd";
	default:
	    return "th";
	}
}		/*end ordinal_suffix*/

/*
*			base_and_dir_name():
*
*	Extracts the base and directory part of a path name.  The
*	strings returned point to substrings of an internal static
*	character strings,  so that the values of these strings are
*	only valid until the next call to this function.
*/

EXPORT	void	base_and_dir_name(
	const char	*name,
	char		**dname,
	char		**bname)
{
	static	char	*buf = NULL;
	static  char    nostring[] = "";
	static	size_t	len = 0;

	if (name == NULL)
	{
	    if (dname != NULL)
	    	*dname = nostring;
	    if (bname != NULL)
	    	*bname = nostring;
	    return;
	}
	if (len == 0)
	{
	    len  = 256 + strlen(name);
	    buf = (char*)malloc(len*sizeof(char));
	}
	else if (strlen(name) >= len)
	{
	    free(buf);
	    len  = 256 + strlen(name);
	    buf = (char*)malloc(len*sizeof(char));
	}
	(void) strcpy(buf,name);
	if (bname != NULL)
	    *bname = basename(buf);
	if (dname != NULL)
	    *dname = dirname(buf);
}		/*end base_and_dir_name*/


EXPORT	void	print_bool(
	const char *mesg1,
	bool	   b,
	const char *mesg2)
{
	if (mesg1 != NULL)
	    (void) printf("%s",mesg1);
	switch (b)
	{
	case YES:
	    (void) printf("YES, FUNCTION_SUCCEEDED");
	    break;
	case NO:
	default:
	    (void) printf("NO, FUNCTION_FAILED");
	    break;
	}
	if (mesg2 != NULL)
	    (void) printf("%s",mesg2);
}		/*end print_bool*/
