/* <<<< Cartographic projection filter program >>>> */
#include "projects.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "emess.h"

/* TK 1999-02-13 */
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__WIN32__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif
/* ! TK 1999-02-13 */

#define MAX_LINE 1000
#define MAX_PARGS 100
#define PJ_INVERS(P) (P->inv ? 1 : 0)
	static PJ
*Proj;
	static projUV
(*proj)(projUV, PJ *);
	static int
reversein = 0,	/* != 0 reverse input arguments */
reverseout = 0,	/* != 0 reverse output arguments */
bin_in = 0,	/* != 0 then binary input */
bin_out = 0,	/* != 0 then binary output */
echoin = 0,	/* echo input data to output line */
tag = '#',	/* beginning of line tag character */
inverse = 0,	/* != 0 then inverse projection */
prescale = 0,	/* != 0 apply cartesian scale factor */
dofactors = 0,	/* determine scale factors */
facs_bad = 0,	/* return condition from pj_factors */
very_verby = 0, /* very verbose mode */
postscale = 0;
	static char
*cheby_str,		/* string controlling Chebychev evaluation */
*oform = (char *)0,	/* output format for x-y or decimal degrees */
*oterr = "*\t*",	/* output line for unprojectable input */
*usage =
"%s\nusage: %s [ -beEfiIlormsStTvVwW [args] ] [ +opts[=arg] ] [ files ]\n";
	static struct FACTORS
facs;
	static double
(*informat)(const char *, char **),	/* input data deformatter function */
fscale = 0.;	/* cartesian scale factor */
	static projUV
int_proj(projUV data) {
	if (prescale) { data.u *= fscale; data.v *= fscale; }
	data = (*proj)(data, Proj);
	if (postscale && data.u != HUGE_VAL)
		{ data.u *= fscale; data.v *= fscale; }
	return(data);
}
	static void	/* file processing function */
process(FILE *fid) {
	char line[MAX_LINE+3], *s, pline[40];
	projUV data;

	for (;;) {
		++emess_dat.File_line;
		if (bin_in) {	/* binary input */
			if (fread(&data, sizeof(projUV), 1, fid) != 1)
				break;
		} else {	/* ascii input */
			if (!(s = fgets(line, MAX_LINE, fid)))
				break;
			if (!strchr(s, '\n')) { /* overlong line */
				int c;
				(void)strcat(s, "\n");
				/* gobble up to newline */
				while ((c = fgetc(fid)) != EOF && c != '\n') ;
			}
			if (*s == tag) {
				if (!bin_out)
					(void)fputs(line, stdout);
				continue;
			}
			if (reversein) {
				data.v = (*informat)(s, &s);
				data.u = (*informat)(s, &s);
			} else {
				data.u = (*informat)(s, &s);
				data.v = (*informat)(s, &s);
			}
			if (data.v == HUGE_VAL)
				data.u = HUGE_VAL;
			if (!*s && (s > line)) --s; /* assumed we gobbled \n */
			if (!bin_out && echoin) {
				int t;
				t = *s;
				*s = '\0';
				(void)fputs(line, stdout);
				*s = t;
				putchar('\t');
			}
		}
		if (data.u != HUGE_VAL) {
			if (prescale) { data.u *= fscale; data.v *= fscale; }
			if (dofactors && !inverse)
				facs_bad = pj_factors(data, Proj, 0., &facs);
			data = (*proj)(data, Proj);
			if (dofactors && inverse)
				facs_bad = pj_factors(data, Proj, 0., &facs);
			if (postscale && data.u != HUGE_VAL)
				{ data.u *= fscale; data.v *= fscale; }
		}
		if (bin_out) { /* binary output */
			(void)fwrite(&data, sizeof(projUV), 1, stdout);
			continue;
		} else if (data.u == HUGE_VAL) /* error output */
			(void)fputs(oterr, stdout);
		else if (inverse && !oform) {	/*ascii DMS output */
			if (reverseout) {
				(void)fputs(rtodms(pline, data.v, 'N', 'S'), stdout);
				putchar('\t');
				(void)fputs(rtodms(pline, data.u, 'E', 'W'), stdout);
			} else {
				(void)fputs(rtodms(pline, data.u, 'E', 'W'), stdout);
				putchar('\t');
				(void)fputs(rtodms(pline, data.v, 'N', 'S'), stdout);
			}
		} else {	/* x-y or decimal degree ascii output */
			if (inverse) {
				data.v *= RAD_TO_DEG;
				data.u *= RAD_TO_DEG;
			}
			if (reverseout) {
				(void)printf(oform,data.v); putchar('\t');
				(void)printf(oform,data.u);
			} else {
				(void)printf(oform,data.u); putchar('\t');
				(void)printf(oform,data.v);
			}
		}
		if (dofactors) /* print scale factor data */
                {
			if (!facs_bad)
				(void)printf("\t<%g %g %g %g %g %g>",
					facs.h, facs.k, facs.s,
					facs.omega * RAD_TO_DEG, facs.a, facs.b);
			else
				(void)fputs("\t<* * * * * *>", stdout);
                }
		(void)fputs(bin_in ? "\n" : s, stdout);
	}
}
	static void	/* file processing function --- verbosely */
vprocess(FILE *fid) {
	char line[MAX_LINE+3], *s, pline[40];
	projUV dat_ll, dat_xy;
	int linvers;

	if (!oform)
		oform = "%.3f";
	if (bin_in || bin_out)
		emess(1,"binary I/O not available in -V option");	
	for (;;) {
		++emess_dat.File_line;
		if (!(s = fgets(line, MAX_LINE, fid)))
			break;
		if (!strchr(s, '\n')) { /* overlong line */
			int c;
			(void)strcat(s, "\n");
			/* gobble up to newline */
			while ((c = fgetc(fid)) != EOF && c != '\n') ;
		}
		if (*s == tag) { /* pass on data */
			(void)fputs(s, stdout);
			continue;
		}
		/* check to override default input mode */
		if (*s == 'I' || *s == 'i') {
			linvers = 1;
			++s;
		} else if (*s == 'I' || *s == 'i') {
			linvers = 0;
			++s;
		} else
			linvers = inverse;
		if (linvers) {
			if (!PJ_INVERS(Proj)) {
				emess(-1,"inverse for this projection not avail.\n");
				continue;
			}
			dat_xy.u = strtod(s, &s);
			dat_xy.v = strtod(s, &s);
			if (dat_xy.u == HUGE_VAL || dat_xy.v == HUGE_VAL) {
				emess(-1,"lon-lat input conversion failure\n");
				continue;
			}
			if (prescale) { dat_xy.u *= fscale; dat_xy.v *= fscale; }
			dat_ll = pj_inv(dat_xy, Proj);
		} else {
			dat_ll.u = dmstor(s, &s);
			dat_ll.v = dmstor(s, &s);
			if (dat_ll.u == HUGE_VAL || dat_ll.v == HUGE_VAL) {
				emess(-1,"lon-lat input conversion failure\n");
				continue;
			}
			dat_xy = pj_fwd(dat_ll, Proj);
			if (postscale) { dat_xy.u *= fscale; dat_xy.v *= fscale; }
		}
		if (pj_errno) {
			emess(-1, pj_strerrno(pj_errno));
			continue;
		}
		if (!*s && (s > line)) --s; /* assumed we gobbled \n */
		if (pj_factors(dat_ll, Proj, 0., &facs)) {
			emess(-1,"failed to conpute factors\n\n");
			continue;
		}
		if (*s != '\n')
			(void)fputs(s, stdout);
		(void)fputs("Longitude: ", stdout);
		(void)fputs(rtodms(pline, dat_ll.u, 'E', 'W'), stdout);
		(void)printf(" [ %.11g ]\n", dat_ll.u * RAD_TO_DEG);
		(void)fputs("Latitude:  ", stdout);
		(void)fputs(rtodms(pline, dat_ll.v, 'N', 'S'), stdout);
		(void)printf(" [ %.11g ]\n", dat_ll.v * RAD_TO_DEG);
		(void)fputs("Easting (x):   ", stdout);
		(void)printf(oform, dat_xy.u); putchar('\n');
		(void)fputs("Northing (y):  ", stdout);
		(void)printf(oform, dat_xy.v); putchar('\n');
		(void)printf("Meridian scale (h)%c: %.8f  ( %.4g %% error )\n",
			facs.code & IS_ANAL_HK ? '*' : ' ', facs.h, (facs.h-1.)*100.);
		(void)printf("Parallel scale (k)%c: %.8f  ( %.4g %% error )\n",
			facs.code & IS_ANAL_HK ? '*' : ' ', facs.k, (facs.k-1.)*100.);
		(void)printf("Areal scale (s):     %.8f  ( %.4g %% error )\n",
			facs.s, (facs.s-1.)*100.);
		(void)printf("Angular distortion (w): %.3f\n", facs.omega *
			RAD_TO_DEG);
		(void)printf("Meridian/Parallel angle: %.5f\n",
			facs.thetap * RAD_TO_DEG);
		(void)printf("Convergence%c: ",facs.code & IS_ANAL_CONV ? '*' : ' ');
		(void)fputs(rtodms(pline, facs.conv, 0, 0), stdout);
		(void)printf(" [ %.8f ]\n", facs.conv * RAD_TO_DEG);
		(void)printf("Max-min (Tissot axis a-b) scale error: %.5f %.5f\n\n",
			facs.a, facs.b);
	}
}

