/* Advance past already accepted output times after a failed/retried increment.
 * CCX uses a one-based pointer into an interleaved (time,value) amplitude array.
 * Only machine-roundoff proximity counts; a genuine future point is retained.
 */
#include <float.h>
#include <math.h>
static int reentrant_advance_completed_points(ITG *next, ITG end,
                                               const double *points,
                                               double now, double period){
    double eps=16*DBL_EPSILON*fmax(fabs(now),fabs(period));
    int count=0;
    while(*next<=end && points[2*(*next)-2]<=now+eps){++*next;++count;}
    return count;
}
