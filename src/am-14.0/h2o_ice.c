/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o_ice.c                  S. Paine rev. 2019 September 26
*
* Water ice properties and absorption.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_types.h"
#include "errlog.h"
#include "h2o_ice.h"
#include "math_const.h"
#include "molecules.h"
#include "phys_const.h"
#include "rayleigh.h"


/*
 * Temperature limits for water ice.
 *
 * A warning threshold is set to inform the user if ice is
 * encountered on a model layer having a temperature exceeding
 * the triple point.  Above this, an error threshold is set to
 * prevent unreasonable extrapolation of density or permittivity
 * models.
 */
static const double H2O_ICE_HIGH_T_ERROR_LIMIT   = 300.;
static const double H2O_ICE_HIGH_T_WARNING_LIMIT = 273.16;


/*
 * The following interpolation table for the complex refractive
 * index of ice at 266 K is adapted from:
 *
 *   S. G. Warren and R. E. Brandt 2008, "Optical constants
 *   of ice from the ultraviolet to the microwave: A revised
 *   compilation."  J. Geophys. Res. 113:D14220.
 *
 * The compilation described therein is available at:
 *
 *   http://www.atmos.washington.edu/ice_optical_constants/
 *
 * The original table gives lambda[um], mr, mi, where 
 *
 *   m = mr + i * mi
 *
 * is the complex refractive index.  The authors state that
 * mr and log(mi) should each be interpolated linearly in
 * log(lambda).
 *
 * For efficiency, the table here has been converted from the
 * original to log(f[GHz]), mr, log(mi), and sorted into
 * ascending frequency order.  Note that linear interpolation in
 * log(f) is equivalent to linear interpolation in log(lambda).
 */
static const struct {
    double log_f;
    double mr;
    double log_mi;
} IOP_tab[] = {
    {  -1.8978,   1.7861,  -7.3239 },
    {  -1.2047,   1.7861,  -8.0020 },
    {  -0.5885,   1.7861,  -8.5711 },
    {  -0.1259,   1.7861,  -8.9388 },
    {   0.2225,   1.7861,  -9.1464 },
    {   0.5101,   1.7861,  -9.2498 },
    {   0.7614,   1.7861,  -9.2765 },
    {   1.0026,   1.7861,  -9.2434 },
    {   1.2487,   1.7861,  -9.1540 },
    {   1.5922,   1.7861,  -8.9526 },
    {   2.0395,   1.7861,  -8.6011 },
    {   2.7587,   1.7861,  -7.9367 },
    {   4.0937,   1.7861,  -6.6173 },
    {   5.4407,   1.7868,  -5.2643 },
    {   6.3962,   1.7908,  -4.2651 },
    {   6.9071,   1.7989,  -3.6695 },
    {   7.2637,   1.8114,  -3.1507 },
    {   7.5357,   1.8268,  -2.7055 },
    {   7.8234,   1.8499,  -2.1620 },
    {   8.0057,   1.8654,  -1.7684 },
    {   8.0552,   1.8681,  -1.6539 },
    {   8.1022,   1.8698,  -1.5418 },
    {   8.1472,   1.8699,  -1.4313 },
    {   8.1682,   1.8688,  -1.3787 },
    {   8.1924,   1.8669,  -1.3171 },
    {   8.2288,   1.8617,  -1.2228 },
    {   8.2763,   1.8504,  -1.0966 },
    {   8.2934,   1.8419,  -1.0376 },
    {   8.3204,   1.8185,  -0.9597 },
    {   8.3548,   1.7780,  -0.8943 },
    {   8.3767,   1.7466,  -0.8754 },
    {   8.3972,   1.7167,  -0.8754 },
    {   8.4129,   1.6945,  -0.8879 },
    {   8.4276,   1.6754,  -0.9103 },
    {   8.4415,   1.6596,  -0.9411 },
    {   8.4554,   1.6471,  -0.9814 },
    {   8.4757,   1.6372,  -1.0581 },
    {   8.4959,   1.6396,  -1.1565 },
    {   8.5000,   1.6419,  -1.1790 },
    {   8.5165,   1.6610,  -1.2783 },
    {   8.5284,   1.6930,  -1.3591 },
    {   8.5435,   1.7410,  -1.3246 },
    {   8.5598,   1.7805,  -1.2103 },
    {   8.5678,   1.7940,  -1.1407 },
    {   8.5784,   1.8041,  -1.0475 },
    {   8.5901,   1.8052,  -0.9519 },
    {   8.6035,   1.7955,  -0.8642 },
    {   8.6174,   1.7792,  -0.8056 },
    {   8.6312,   1.7687,  -0.7767 },
    {   8.6453,   1.7668,  -0.7381 },
    {   8.6735,   1.7649,  -0.6249 },
    {   8.6982,   1.7505,  -0.4795 },
    {   8.7235,   1.6841,  -0.2992 },
    {   8.7396,   1.5874,  -0.2036 },
    {   8.7531,   1.4725,  -0.1675 },
    {   8.7669,   1.3474,  -0.1919 },
    {   8.7789,   1.2543,  -0.2634 },
    {   8.7931,   1.1853,  -0.3883 },
    {   8.8136,   1.1496,  -0.5915 },
    {   8.8362,   1.1549,  -0.7939 },
    {   8.8552,   1.1742,  -0.9301 },
    {   8.8939,   1.2054,  -1.1270 },
    {   8.9320,   1.2243,  -1.2542 },
    {   8.9670,   1.2317,  -1.3752 },
    {   8.9999,   1.2371,  -1.5019 },
    {   9.0374,   1.2445,  -1.6634 },
    {   9.0705,   1.2530,  -1.8233 },
    {   9.1044,   1.2641,  -2.0040 },
    {   9.1769,   1.2969,  -2.4517 },
    {   9.2585,   1.3404,  -2.9716 },
    {   9.3527,   1.3854,  -3.3814 },
    {   9.3920,   1.4030,  -3.5066 },
    {   9.4081,   1.4102,  -3.5439 },
    {   9.4496,   1.4282,  -3.6009 },
    {   9.4732,   1.4390,  -3.6119 },
    {   9.4929,   1.4486,  -3.5936 },
    {   9.5099,   1.4575,  -3.5405 },
    {   9.5745,   1.4877,  -3.1011 },
    {   9.5948,   1.4956,  -2.9004 },
    {   9.6151,   1.4986,  -2.7031 },
    {   9.6348,   1.4963,  -2.5903 },
    {   9.6544,   1.4939,  -2.5770 },
    {   9.6733,   1.4982,  -2.6311 },
    {   9.6871,   1.5043,  -2.5963 },
    {   9.7029,   1.5104,  -2.5572 },
    {   9.7105,   1.5141,  -2.5383 },
    {   9.7636,   1.5303,  -2.2349 },
    {   9.7806,   1.5306,  -2.1542 },
    {   9.7972,   1.5294,  -2.0956 },
    {   9.8142,   1.5302,  -2.0794 },
    {   9.8302,   1.5353,  -2.0557 },
    {   9.8617,   1.5481,  -1.9519 },
    {   9.8778,   1.5559,  -1.8839 },
    {   9.8928,   1.5637,  -1.8079 },
    {   9.9223,   1.5762,  -1.6195 },
    {   9.9513,   1.5775,  -1.4024 },
    {   9.9661,   1.5701,  -1.3056 },
    {   9.9797,   1.5596,  -1.2242 },
    {   9.9934,   1.5458,  -1.1552 },
    {  10.0074,   1.5300,  -1.0936 },
    {  10.0208,   1.5132,  -1.0385 },
    {  10.0337,   1.4928,  -0.9835 },
    {  10.0467,   1.4683,  -0.9442 },
    {  10.0598,   1.4412,  -0.9088 },
    {  10.0851,   1.3822,  -0.8627 },
    {  10.1094,   1.3194,  -0.8627 },
    {  10.1343,   1.2546,  -0.8940 },
    {  10.1573,   1.1983,  -0.9702 },
    {  10.1807,   1.1439,  -1.0759 },
    {  10.2030,   1.1023,  -1.2730 },
    {  10.2130,   1.0886,  -1.3943 },
    {  10.2248,   1.0833,  -1.5896 },
    {  10.2359,   1.0867,  -1.7838 },
    {  10.2462,   1.0971,  -2.0099 },
    {  10.2566,   1.1136,  -2.2256 },
    {  10.2671,   1.1323,  -2.4304 },
    {  10.2777,   1.1501,  -2.5903 },
    {  10.2885,   1.1659,  -2.7394 },
    {  10.3083,   1.1926,  -2.9941 },
    {  10.3281,   1.2151,  -3.1123 },
    {  10.3570,   1.2387,  -3.2124 },
    {  10.3852,   1.2561,  -3.2731 },
    {  10.4036,   1.2655,  -3.2952 },
    {  10.4215,   1.2735,  -3.3066 },
    {  10.4480,   1.2835,  -3.3093 },
    {  10.4737,   1.2917,  -3.3031 },
    {  10.4906,   1.2964,  -3.2968 },
    {  10.5071,   1.3007,  -3.2904 },
    {  10.5233,   1.3047,  -3.2808 },
    {  10.5393,   1.3086,  -3.2652 },
    {  10.5552,   1.3123,  -3.2414 },
    {  10.5707,   1.3158,  -3.2087 },
    {  10.5859,   1.3183,  -3.1489 },
    {  10.6009,   1.3200,  -3.1145 },
    {  10.6157,   1.3218,  -3.0791 },
    {  10.6304,   1.3229,  -3.0283 },
    {  10.6447,   1.3236,  -3.0078 },
    {  10.6590,   1.3244,  -2.9640 },
    {  10.6798,   1.3242,  -2.9136 },
    {  10.7003,   1.3229,  -2.8682 },
    {  10.7137,   1.3219,  -2.8647 },
    {  10.7270,   1.3214,  -2.8508 },
    {  10.7400,   1.3207,  -2.8304 },
    {  10.7465,   1.3201,  -2.8194 },
    {  10.7530,   1.3194,  -2.8083 },
    {  10.7594,   1.3187,  -2.7967 },
    {  10.7657,   1.3178,  -2.7855 },
    {  10.7720,   1.3167,  -2.7702 },
    {  10.7783,   1.3152,  -2.7551 },
    {  10.7968,   1.3083,  -2.7242 },
    {  10.8090,   1.3015,  -2.7428 },
    {  10.8389,   1.2915,  -2.9868 },
    {  10.8447,   1.2904,  -3.0698 },
    {  10.8506,   1.2902,  -3.1536 },
    {  10.8564,   1.2906,  -3.2486 },
    {  10.8622,   1.2917,  -3.3413 },
    {  10.8679,   1.2933,  -3.4247 },
    {  10.8960,   1.3013,  -3.8023 },
    {  10.9244,   1.3100,  -4.0846 },
    {  10.9501,   1.3176,  -4.2831 },
    {  10.9622,   1.3212,  -4.3413 },
    {  10.9816,   1.3268,  -4.4063 },
    {  11.0014,   1.3325,  -4.3901 },
    {  11.0208,   1.3379,  -4.3073 },
    {  11.0592,   1.3470,  -3.9425 },
    {  11.0891,   1.3491,  -3.5463 },
    {  11.0935,   1.3482,  -3.5042 },
    {  11.0964,   1.3473,  -3.4770 },
    {  11.1057,   1.3444,  -3.4486 },
    {  11.1193,   1.3412,  -3.5099 },
    {  11.1322,   1.3401,  -3.6283 },
    {  11.1411,   1.3406,  -3.7251 },
    {  11.1665,   1.3447,  -3.9808 },
    {  11.2001,   1.3526,  -4.2192 },
    {  11.2323,   1.3623,  -4.4990 },
    {  11.2636,   1.3750,  -4.7978 },
    {  11.2824,   1.3850,  -4.9667 },
    {  11.2939,   1.3924,  -5.0098 },
    {  11.3233,   1.4146,  -4.8102 },
    {  11.3414,   1.4328,  -4.4830 },
    {  11.3484,   1.4411,  -4.3230 },
    {  11.3555,   1.4502,  -4.1440 },
    {  11.3627,   1.4604,  -3.9523 },
    {  11.3696,   1.4710,  -3.7520 },
    {  11.3765,   1.4827,  -3.5450 },
    {  11.3833,   1.4949,  -3.3459 },
    {  11.3900,   1.5086,  -3.1498 },
    {  11.3968,   1.5241,  -2.9352 },
    {  11.4034,   1.5396,  -2.7061 },
    {  11.4103,   1.5559,  -2.4817 },
    {  11.4169,   1.5714,  -2.2730 },
    {  11.4233,   1.5905,  -2.1120 },
    {  11.4300,   1.6108,  -1.8452 },
    {  11.4365,   1.6248,  -1.6451 },
    {  11.4430,   1.6405,  -1.4355 },
    {  11.4493,   1.6477,  -1.2040 },
    {  11.4555,   1.6336,  -0.9889 },
    {  11.4619,   1.5970,  -0.8255 },
    {  11.4682,   1.5478,  -0.7113 },
    {  11.4746,   1.4933,  -0.6180 },
    {  11.4807,   1.4225,  -0.5226 },
    {  11.4869,   1.3215,  -0.4700 },
    {  11.4931,   1.2089,  -0.4910 },
    {  11.4990,   1.1259,  -0.5906 },
    {  11.5053,   1.0722,  -0.7072 },
    {  11.5112,   1.0390,  -0.8255 },
    {  11.5172,   1.0180,  -0.9467 },
    {  11.5233,   1.0026,  -1.0584 },
    {  11.5290,   0.9873,  -1.1648 },
    {  11.5348,   0.9678,  -1.2874 },
    {  11.5410,   0.9538,  -1.5096 },
    {  11.5465,   0.9563,  -1.7779 },
    {  11.5524,   0.9747,  -2.1219 },
    {  11.5583,   1.0001,  -2.4595 },
    {  11.5639,   1.0236,  -2.7667 },
    {  11.5695,   1.0453,  -3.0852 },
    {  11.5752,   1.0657,  -3.4281 },
    {  11.5891,   1.1083,  -4.3080 },
    {  11.6084,   1.1507,  -5.6527 },
    {  11.6269,   1.1776,  -6.5951 },
    {  11.6477,   1.1983,  -7.1236 },
    {  11.6553,   1.2043,  -7.1897 },
    {  11.6592,   1.2071,  -7.2039 },
    {  11.6631,   1.2097,  -7.2099 },
    {  11.6689,   1.2135,  -7.1990 },
    {  11.6748,   1.2169,  -7.1760 },
    {  11.6866,   1.2232,  -7.1711 },
    {  11.6946,   1.2270,  -7.1914 },
    {  11.7107,   1.2337,  -7.2905 },
    {  11.7230,   1.2382,  -7.3763 },
    {  11.7312,   1.2409,  -7.4259 },
    {  11.7396,   1.2435,  -7.4867 },
    {  11.7480,   1.2459,  -7.5642 },
    {  11.7564,   1.2482,  -7.6862 },
    {  11.7650,   1.2504,  -7.8602 },
    {  11.7736,   1.2525,  -8.0668 },
    {  11.7823,   1.2545,  -8.2790 },
    {  11.7911,   1.2564,  -8.4352 },
    {  11.7955,   1.2573,  -8.4789 },
    {  11.7999,   1.2582,  -8.4998 },
    {  11.8021,   1.2587,  -8.5043 },
    {  11.8044,   1.2591,  -8.4989 },
    {  11.8133,   1.2609,  -8.4092 },
    {  11.8269,   1.2633,  -8.2145 },
    {  11.8361,   1.2648,  -8.0363 },
    {  11.8454,   1.2663,  -7.8178 },
    {  11.8547,   1.2677,  -7.5512 },
    {  11.8665,   1.2694,  -7.1840 },
    {  11.8770,   1.2707,  -6.8850 },
    {  11.8872,   1.2718,  -6.6711 },
    {  11.8974,   1.2728,  -6.5307 },
    {  11.9077,   1.2736,  -6.4592 },
    {  11.9177,   1.2744,  -6.4131 },
    {  11.9277,   1.2750,  -6.4290 },
    {  11.9374,   1.2756,  -6.5417 },
    {  11.9471,   1.2762,  -6.8052 },
    {  11.9570,   1.2771,  -7.2559 },
    {  11.9664,   1.2780,  -7.8040 },
    {  11.9743,   1.2788,  -8.2680 },
    {  11.9849,   1.2797,  -8.7205 },
    {  11.9903,   1.2802,  -8.8513 },
    {  11.9930,   1.2805,  -8.8904 },
    {  11.9957,   1.2807,  -8.9184 },
    {  12.0011,   1.2811,  -8.9403 },
    {  12.0065,   1.2816,  -8.9464 },
    {  12.0231,   1.2828,  -8.8660 },
    {  12.0455,   1.2843,  -8.7903 },
    {  12.0627,   1.2853,  -8.7096 },
    {  12.0802,   1.2863,  -8.5817 },
    {  12.0921,   1.2870,  -8.4945 },
    {  12.1101,   1.2879,  -8.3513 },
    {  12.1327,   1.2890,  -8.2324 },
    {  12.1490,   1.2897,  -8.0773 },
    {  12.1642,   1.2903,  -7.8602 },
    {  12.1804,   1.2909,  -7.6856 },
    {  12.1875,   1.2912,  -7.6195 },
    {  12.1954,   1.2914,  -7.5727 },
    {  12.2027,   1.2916,  -7.5290 },
    {  12.2101,   1.2918,  -7.4998 },
    {  12.2181,   1.2920,  -7.6128 },
    {  12.2249,   1.2921,  -7.8273 },
    {  12.2324,   1.2924,  -8.1313 },
    {  12.2400,   1.2927,  -8.5023 },
    {  12.2462,   1.2929,  -8.7943 },
    {  12.2532,   1.2931,  -9.1827 },
    {  12.2602,   1.2934,  -9.7280 },
    {  12.2673,   1.2937, -10.2769 },
    {  12.2744,   1.2939, -10.8298 },
    {  12.2815,   1.2941, -10.9590 },
    {  12.2888,   1.2944, -11.0555 },
    {  12.2960,   1.2946, -11.1209 },
    {  12.3034,   1.2949, -11.1623 },
    {  12.3107,   1.2951, -11.1836 },
    {  12.3182,   1.2953, -11.2203 },
    {  12.3257,   1.2955, -11.2353 },
    {  12.3332,   1.2957, -11.2353 },
    {  12.3408,   1.2959, -11.2429 },
    {  12.3485,   1.2961, -11.2353 },
    {  12.3562,   1.2963, -11.2353 },
    {  12.3640,   1.2965, -11.2277 },
    {  12.3718,   1.2967, -11.2128 },
    {  12.3797,   1.2969, -11.2353 },
    {  12.3877,   1.2971, -11.2583 },
    {  12.3957,   1.2973, -11.3141 },
    {  12.4038,   1.2975, -11.3907 },
    {  12.4120,   1.2977, -11.4931 },
    {  12.4202,   1.2979, -11.6568 },
    {  12.4285,   1.2980, -11.9119 },
    {  12.4369,   1.2982, -12.0646 },
    {  12.4453,   1.2984, -12.2532 },
    {  12.4538,   1.2986, -12.4700 },
    {  12.4624,   1.2988, -12.7037 },
    {  12.4711,   1.2990, -12.9870 },
    {  12.4798,   1.2991, -13.0046 },
    {  12.4886,   1.2993, -13.1026 },
    {  12.4975,   1.2995, -13.2167 },
    {  12.5065,   1.2997, -13.2502 },
    {  12.5155,   1.2998, -13.2849 },
    {  12.5247,   1.3000, -13.2674 },
    {  12.5339,   1.3002, -13.2616 },
    {  12.5432,   1.3003, -13.2222 },
    {  12.5526,   1.3005, -13.1426 },
    {  12.5621,   1.3007, -13.0408 },
    {  12.5716,   1.3009, -12.9696 },
    {  12.5813,   1.3010, -12.9696 },
    {  12.5910,   1.3012, -13.0046 },
    {  12.6009,   1.3014, -13.1224 },
    {  12.6108,   1.3015, -13.3331 },
    {  12.6209,   1.3017, -13.5303 },
    {  12.6310,   1.3019, -13.7022 },
    {  12.6413,   1.3020, -13.8924 },
    {  12.6517,   1.3022, -14.0965 },
    {  12.6621,   1.3023, -14.3230 },
    {  12.6727,   1.3025, -14.4079 },
    {  12.6834,   1.3027, -14.4869 },
    {  12.6942,   1.3028, -14.5621 },
    {  12.7052,   1.3030, -14.6274 },
    {  12.7162,   1.3032, -14.6830 },
    {  12.7274,   1.3033, -14.7520 },
    {  12.7387,   1.3035, -14.9091 },
    {  12.7501,   1.3037, -15.1435 },
    {  12.7617,   1.3039, -15.3526 },
    {  12.7734,   1.3040, -15.5138 },
    {  12.7852,   1.3042, -15.7060 },
    {  12.7972,   1.3044, -15.7465 },
    {  12.8093,   1.3046, -15.7604 },
    {  12.8216,   1.3047, -15.7816 },
    {  12.8340,   1.3049, -15.8254 },
    {  12.8466,   1.3051, -15.9526 },
    {  12.8593,   1.3053, -16.0983 },
    {  12.8722,   1.3055, -16.2712 },
    {  12.8853,   1.3057, -16.4634 },
    {  12.8985,   1.3059, -16.6508 },
    {  12.9120,   1.3060, -16.8274 },
    {  12.9256,   1.3062, -16.9621 },
    {  12.9393,   1.3065, -17.0269 },
    {  12.9533,   1.3067, -17.1852 },
    {  12.9675,   1.3069, -17.3560 },
    {  12.9819,   1.3071, -17.5452 },
    {  12.9965,   1.3073, -17.6835 },
    {  13.0113,   1.3076, -17.7841 },
    {  13.0264,   1.3078, -17.9139 },
    {  13.0416,   1.3080, -18.0630 },
    {  13.0571,   1.3083, -18.2218 },
    {  13.0729,   1.3085, -18.3815 },
    {  13.0889,   1.3088, -18.5738 },
    {  13.1051,   1.3091, -18.7932 },
    {  13.1217,   1.3094, -18.9776 },
    {  13.1385,   1.3097, -19.1279 },
    {  13.1556,   1.3100, -19.2980 },
    {  13.1730,   1.3103, -19.4817 },
    {  13.1907,   1.3106, -19.6798 },
    {  13.2087,   1.3110, -19.8952 },
    {  13.2270,   1.3114, -20.1283 },
    {  13.2457,   1.3117, -20.3804 },
    {  13.2648,   1.3121, -20.6500 },
    {  13.2842,   1.3126, -20.9419 },
    {  13.3040,   1.3130, -21.2528 },
    {  13.3242,   1.3135, -21.5975 },
    {  13.3448,   1.3140, -21.9747 },
    {  13.3659,   1.3145, -22.3549 },
    {  13.3874,   1.3151, -22.7444 },
    {  13.4094,   1.3157, -23.1050 },
    {  13.4318,   1.3163, -23.4930 },
    {  13.4548,   1.3170, -23.9077 },
    {  13.4783,   1.3177, -24.1858 },
    {  13.5024,   1.3185, -24.3467 },
    {  13.5271,   1.3194, -24.4677 },
    {  13.5525,   1.3203, -24.6353 },
    {  13.6607,   1.3249, -24.6353 },
    {  13.8148,   1.3339, -24.6353 },
    {  13.9971,   1.3509, -24.6353 },
    {  14.1715,   1.3801, -24.6353 },
    {  14.2108,   1.3901, -24.6353 },
    {  14.2153,   1.3914, -24.1501 },
    {  14.2253,   1.3943, -23.0703 },
    {  14.2354,   1.3974, -21.9905 },
    {  14.2456,   1.4007, -20.9109 },
    {  14.2559,   1.4043, -19.8313 },
    {  14.2663,   1.4081, -18.7518 },
    {  14.2769,   1.4122, -17.6726 },
    {  14.2875,   1.4167, -16.5932 },
    {  14.2982,   1.4215, -15.5143 },
    {  14.3091,   1.4268, -14.4311 },
    {  14.3201,   1.4326, -13.3644 },
    {  14.3312,   1.4390, -12.4092 },
    {  14.3425,   1.4462, -11.4537 },
    {  14.3538,   1.4543, -10.4988 },
    {  14.3653,   1.4635,  -9.5437 },
    {  14.3769,   1.4743,  -8.5887 },
    {  14.3887,   1.4870,  -7.6340 },
    {  14.4006,   1.5027,  -6.6798 },
    {  14.4127,   1.5225,  -5.7251 },
    {  14.4249,   1.5482,  -4.7709 },
    {  14.4372,   1.5818,  -3.8167 },
    {  14.4503,   1.6214,  -2.8473 },
    {  14.4560,   1.6309,  -2.4889 },
    {  14.4630,   1.6353,  -2.2073 },
    {  14.4752,   1.6363,  -1.8579 },
    {  14.4816,   1.6345,  -1.7603 },
    {  14.4875,   1.6311,  -1.6195 },
    {  14.5000,   1.6220,  -1.4784 },
    {  14.5120,   1.6156,  -1.2946 },
    {  14.5180,   1.6077,  -1.2107 },
    {  14.5214,   1.6021,  -1.1552 },
    {  14.5241,   1.5915,  -1.1026 },
    {  14.5295,   1.5729,  -1.0642 },
    {  14.5357,   1.5484,  -0.9755 },
    {  14.5419,   1.5157,  -0.9442 },
    {  14.5474,   1.4845,  -0.9039 },
    {  14.5523,   1.4525,  -0.9002 },
    {  14.5593,   1.4130,  -0.9138 },
    {  14.5706,   1.3575,  -0.9623 },
    {  14.5820,   1.3151,  -1.0671 },
    {  14.5863,   1.2970,  -1.0966 },
    {  14.5928,   1.2851,  -1.2413 },
    {  14.5986,   1.2849,  -1.3056 },
    {  14.6045,   1.2871,  -1.4147 },
    {  14.6148,   1.3030,  -1.5465 },
    {  14.6260,   1.3210,  -1.6399 },
    {  14.6366,   1.3389,  -1.7148 },
    {  14.6472,   1.3607,  -1.7603 },
    {  14.6549,   1.3777,  -1.7487 },
    {  14.6728,   1.3999,  -1.6348 },
    {  14.6983,   1.4061,  -1.4961 },
    {  14.7228,   1.3929,  -1.3984 },
    {  14.7471,   1.3783,  -1.4106 },
    {  14.7824,   1.3772,  -1.4313 },
    {  14.7939,   1.3806,  -1.4313 },
    {  14.8181,   1.3886,  -1.4106 },
    {  14.8383,   1.3981,  -1.3823 },
    {  14.8810,   1.4047,  -1.2344 },
    {  14.9215,   1.3964,  -1.1457 },
    {  14.9605,   1.3948,  -1.0818 },
    {  14.9860,   1.4028,  -1.0217 },
    {  14.9990,   1.4091,  -0.9676 },
    {  15.0345,   1.3807,  -0.7593 },
    {  15.0701,   1.3117,  -0.6951 },
    {  15.1034,   1.2629,  -0.6931 },
    {  15.1366,   1.2273,  -0.6733 },
    {  15.1683,   1.1885,  -0.6463 },
    {  15.1998,   1.1449,  -0.6330 },
    {  15.2172,   1.1183,  -0.6274 },
    {  15.2295,   1.0986,  -0.6292 },
    {  15.2587,   1.0536,  -0.6387 },
    {  15.2859,   1.0093,  -0.6597 },
    {  15.3139,   0.9679,  -0.7093 },
    {  15.3396,   0.9400,  -0.7679 },
    {  15.3660,   0.9173,  -0.8210 },
    {  15.3915,   0.8970,  -0.8771 },
    {  15.4160,   0.8797,  -0.9365 },
    {  15.4411,   0.8647,  -1.0051 },
    {  15.4633,   0.8566,  -1.0788 },
    {  15.4861,   0.8519,  -1.1426 },
    {  15.5095,   0.8489,  -1.2140 },
    {  15.5315,   0.8497,  -1.2765 },
    {  15.5521,   0.8505,  -1.3168 },
    {  15.5731,   0.8483,  -1.3509 },
    {  15.5946,   0.8428,  -1.3863 },
    {  15.6146,   0.8347,  -1.4271 },
    {  15.6350,   0.8281,  -1.5006 },
    {  15.6537,   0.8263,  -1.5702 },
    {  15.6727,   0.8258,  -1.6348 },
    {  15.6921,   0.8255,  -1.6983 },
    {  15.7097,   0.8250,  -1.7545 },
    {  15.7276,   0.8228,  -1.8079 },
};

#define IOP_NTAB    (sizeof(IOP_tab) / sizeof(IOP_tab[0]))

double H2O_ice_density(double);
int H2O_ice_permittivity(
        double*, double*, const double*, const gridsize_t, const double);
int H2O_ice_permittivity_Maetzler(
        double*, double*, const double*, const gridsize_t, const double);
int H2O_ice_permittivity_Warren_Brandt(
        double*, double*, const double*, const gridsize_t);


/***********************************************************
* int iwp_abs_Rayleigh(
*         double *k,
*         const double *f,
*         const gridsize_t ngrid,
*         const double T)
*
* Purpose:
*   Computes the molecular absorption coefficient [cm^2] for
*   water ice particles, in the Rayleigh limit.
*
*   Internally, ice water path (IWP) is represented as a
*   molecular column density in [cm^-2], but configuration
*   files will typically specify IWP in practical units such
*   as [kg*m^-2] or [mm_pwv].
*
* Arguments:
*   double *k - absorption coefficient [cm^2]
*   const double *f - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int iwp_abs_Rayleigh(
        double *k,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    double *epsr = NULL;
    double *epsi = NULL;
    double rho;

    /*
     * Check temperature limits.
     */
    if (T > H2O_ICE_HIGH_T_ERROR_LIMIT) {
        errlog(158, 0);
        return 1;
    }
    if (T > H2O_ICE_HIGH_T_WARNING_LIMIT) {
        errlog(159, 0);
    }
    /*
     * Allocate temporary arrays and compute the complex permittivity
     * at temperature T.
     */
    if ((epsr = (double*)malloc(ngrid * sizeof(double))) == NULL) {
        errlog(143, 0);
        return 1;
    }
    if ((epsi = (double*)malloc(ngrid * sizeof(double))) == NULL) {
        free(epsr);
        errlog(143, 0);
        return 1;
    }
    if (H2O_ice_permittivity(epsr, epsi, f, ngrid, T)) {
        free(epsi);
        free(epsr);
        return 1;
    }
    /*
     * H2O_ice_density() returns [molecules/cm^3].
     */
    rho = H2O_ice_density(T);
    Rayleigh_mol_abscoeff(k, f, ngrid, epsr, epsi, rho);
    free(epsi);
    free(epsr);
    return 0;
}   /* iwp_abs_Rayleigh() */


/***********************************************************
* double H2O_ice_density(double T)
*
* Purpose:
*   Computes the density of water ice [molecules / cm^3],
*   as a function of temperature at atmospheric pressure,
*   using the fit to x-ray diffraction measurements in
*
*     K. Roettger, A. Endriss, J. Ihringer, S. Doyle, and
*     W. F. Kuhs 1994, "Lattice Constants and Thermal
*     Expansion of H2O and D2O Ice Ih Between 10 and 265K."
*     Acta Cryst. B50:644
*
*   and the Addendum by the same authors, giving higher-
*   precision coefficients:
*
*     K. Roettger, A. Endriss, J. Ihringer, S. Doyle, and
*     W. F. Kuhs 2012, "Lattice Constants and Thermal
*     Expansion of H2O and D2O Ice Ih Between 10 and 265K.
*     Addendum."  Acta Cryst. B68:91
*
* Arguments:
*   double T - Temperature [K], between 10 K and 273 K.
*     This argument is assumed to have been checked by the
*     calling function.
*
* Return:
*   Density [molecules / cm^3].
************************************************************/

double H2O_ice_density(double T)
{
    static const double a0 = 128.2147;
    static const double a3 = -1.3152e-6;
    static const double a4 =  2.4837e-8;
    static const double a5 = -1.6064e-10;
    static const double a6 =  4.6097e-13;
    static const double a7 = -4.9661e-16;
    
    /*
     * Here, v is the unit cell volume in angstrom^3, where
     * a unit cell of ice 1h comprises four molecules.  The
     * density in molecules / cm3 is thus 4 * 1e24 / v.
     */
    double v = a0 + T * T * T * (a3 + T * (a4 + T * (a5 + T * (a6 + T * a7))));
    return 4.0e24 / v;
}   /* H2O_ice_density() */


/***********************************************************
* int H2O_ice_permittivity(
*         double *epsr,
*         double *epsi,
*         const double *f,
*         const gridsize_t ngrid,
*         const double T)
*
* Purpose:
*   Computes the complex permittivity of ice, by splicing
*   together H2O_ice_permittivity_Maetzler() at microwave
*   frequencies, and H2O_ice_permittivity_Warren_Brandt()
*   at millimeter-wave to optical frequencies.
*
*   Although this function will accept frequencies down to
*   f = 0, it will not give correct results near and below
*   the lowest relaxation frequency in ice, which ranges from
*   about 100 kHz downwards, depending on temperature.  A
*   good reference is
*
*     V. F. Petrenko and R. W. Whitworth 2002 "Physics of
*     Ice." Oxford University Press.
*   
* Arguments:
*   double *epsr - pointer to array to receive the real part
*                  of the permittivity
*   double *epsi - pointer to array to receive the imaginary
*                  part of the permittivity
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T         - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int H2O_ice_permittivity(
        double *epsr,
        double *epsi,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    /*
     * The temperature-dependent Maetzler permittivity and the
     * temperature-independent Warren-Brandt permittivity are
     * spliced together at frequency fsplice, defined here.  This
     * is done by adjusting the Warren-Brandt values, applying a
     * shift to the real part and a scaling to the imaginary
     * part.
     *
     * Besides maintaining continuity, this adjustment serves as
     * an estimate of the temperature dependence for f > fsplice.
     * See the discussion in
     *
     *   K. Woschnagg and P. B. Price 2001, "Temperature
     *   dependence of absorption in ice at 532 nm." Applied
     *   Optics 40:2496
     *
     * on the approximate scaling behavior of ice absorption from
     * microwave through optical frequencies.
     */
    static const double fsplice = 100.; /* GHz */
    gridsize_t m; /* number of grid points below fsplice */

    /*
     * To find m, assume the frequency grid is in ascending order,
     * though not necessarily uniform.
     */
    if (fsplice < f[0]) {
        m = 0;
    } else if (fsplice > f[ngrid - 1]) {
        m = ngrid;
    } else {
        gridsize_t ilow = 0;
        gridsize_t imid;
        m = ngrid - 1;
        while (m - ilow > 1) {
            imid = (ilow + m) /2;
            if (fsplice >= f[imid])
                ilow = imid;
            else
                m = imid;
        }
    }
    if (m > 0) {
        if (H2O_ice_permittivity_Maetzler(epsr, epsi, f, ngrid, T))
            return 1;
    }
    if (m < ngrid) {
        double epsr_M, epsi_M, epsr_WB, epsi_WB;
        double epsr_shift, epsi_scale;
        if (H2O_ice_permittivity_Maetzler(&epsr_M, &epsi_M, &fsplice, 1, T))
            return 1;
        if (H2O_ice_permittivity_Warren_Brandt(
            &epsr_WB, &epsi_WB, &fsplice, 1))
            return 1;
        epsr_shift = epsr_M - epsr_WB;
        epsi_scale = epsi_M / epsi_WB;
        if (H2O_ice_permittivity_Warren_Brandt(
            epsr + m, epsi + m, f + m, ngrid - m))
            return 1;
        for (; m < ngrid; ++m) {
            epsr[m] += epsr_shift;
            epsi[m] *= epsi_scale;
        }
    }
    return 0;
}   /* H2O_ice_permittivity() */


/***********************************************************
* int H2O_ice_permittivity_Maetzler(
*         double *epsr,
*         double *epsi,
*         const double *f,
*         const gridsize_t ngrid,
*         const double T)
*
* Purpose:
*   Computes the complex permittivity of pure water ice
*   at RF through millimeter-wave frequencies, using the
*   formulation in
*
*     C. Maetzler in C. Maetzler, ed. 2006 "Thermal Microwave
*     Radiation: Applications for Remote Sensing." Section 5.3,
*     Institution of Engineering and Technology, London.
*   
* Arguments:
*   double *epsr - pointer to array to receive the real part
*                  of the permittivity
*   double *epsi - pointer to array to receive the imaginary
*                  part of the permittivity
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T         - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int H2O_ice_permittivity_Maetzler(
        double *epsr,
        double *epsi,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    gridsize_t i;
    double er, alpha, beta_M0, dbeta, theta;
    double r0;
    /*
     * Constants for the real part of the permittivity, from
     *
     *   C. Maetzler and U. Wegmueller 1987, "Dielectric
     *   properties of fresh-water ice at microwave frequencies."
     *   J. Phys D: Appl. Phys. 20:1623.
     */
    static const double Tmin_MW  = 243.0;
    static const double er0 = 3.1884;
    static const double er1 = 9.1e-4;
    /*
     * Low-temperature limit of the microwave real permittivity,
     * from
     *
     *   S. R. Gough 1972, "A Low Temperature Dielectric Cell and
     *   the Permittivity of Hexagonal Ice to 2 K."  Canadian J.
     *   of Chemistry 50:3046.
     */
    static const double er_Gough = 3.093;
    /*
     * Constants for the imaginary part of the permittivity, from
     *
     *   C. Maetzler in C. Maetzler, ed. 2006 "Thermal Microwave
     *   Radiation: Applications for Remote Sensing." Section
     *   5.3.3.1, Institution of Engineering and Technology,
     *   London.
     */
    static const double T0     = 300.;
    static const double alpha0 = 0.00504;
    static const double alpha1 = 0.0062;
    static const double c      = -22.1;
    static const double B1     = 0.0207;
    static const double B2     = 1.16e-11;
    static const double T_b    = 335.;
    static const double dbeta0 = -9.963;
    static const double dbeta1 = 0.0372;
    static const double T_triple_point = 273.16;

    /*
     * For the real part of the permittivity, use the Maetzler-
     * Wegmueller linear formula down to 243 K.  At lower
     * temperature, use a Hermite spline between the M-W formula
     * at 243 K and the low T limit of Gough, assuming d(epsr)/dT
     * = 0 at T = 0.
     */
    if (T >= Tmin_MW) {
        er = er0 + er1 * (T - T_STP);
    } else {
        double p = T / Tmin_MW;
        double a0 = p - 1.0;
        double a1 = p * p;
        double b1 = a1;
        a1 *= (3.0 - 2.0 * p);
        b1 *= a0;
        a0 *= a0 * (1.0 + 2.0 * p);
        er = a0 * er_Gough +
             a1 * (er0 + er1 * (Tmin_MW - T_STP)) +
             b1 * (er1 * Tmin_MW);
    }
    for (i = 0; i < ngrid; ++i)
        epsr[i] = er;
    /*
     * Imaginary part.
     */
    if (T0 > T * (1.0 + (DBL_MIN_10_EXP * LOG_10) / c)) {
        alpha = 0.0;
    } else {
        theta = (T0 / T) - 1.0;
        alpha = (alpha0 + alpha1 * theta) * exp(c * theta);
    }
    if (T_b > T * DBL_MAX_10_EXP * LOG_10) {
        beta_M0 = 0.0;
    } else {
        r0 = exp(T_b / T);
        beta_M0 = B1 * r0;
        r0 -= 1.0;
        r0 *= r0;
        beta_M0 /= T * r0;
    }
    dbeta = exp(dbeta0 + dbeta1 * (T - T_triple_point));
    for (i = 0; i < ngrid; ++i) {
        double rf = f[i] + F_EPSILON;
        double beta_M = beta_M0 + B2 * rf * rf; 
        double beta = beta_M + dbeta;
        epsi[i] = alpha / rf + beta * rf;
    }
    return 0;
}   /* H2O_ice_permittivity_Maetzler() */


/***********************************************************
* int H2O_ice_permittivity_Warren_Brandt(
*         double *epsr,
*         double *epsi,
*         const double *f,
*         const gridsize_t ngrid)
*
* Purpose:
*   Computes the complex permittivity of pure water ice,
*   using the tabulated data of
*
*     S. G. Warren and R. E. Brandt 2008, "Optical constants
*     of ice from the ultraviolet to the microwave: A revised
*     compilation."  J. Geophys. Res. 113:D14220.
*   
* Arguments:
*   double *epsr - pointer to array to receive the real part
*                  of the permittivity
*   double *epsi - pointer to array to receive the imaginary
*                  part of the permittivity
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int H2O_ice_permittivity_Warren_Brandt(
        double *epsr,
        double *epsi,
        const double *f,
        const gridsize_t ngrid)
{
    gridsize_t i;
    unsigned int j = 0;
    double log_f1  = IOP_tab[j].log_f;
    double log_f2  = IOP_tab[j+1].log_f;
    double dlog_f  = log_f2 - log_f1;
    double mr1     = IOP_tab[j].mr;
    double mr2     = IOP_tab[j+1].mr;
    double log_mi1 = IOP_tab[j].log_mi;
    double log_mi2 = IOP_tab[j+1].log_mi;

    if (log(f[0] + F_EPSILON) < IOP_tab[0].log_f ||
        log(f[ngrid - 1] + F_EPSILON) > IOP_tab[IOP_NTAB - 1].log_f) {
        errlog(142, 0);
        return 1;
    }
    for (i = 0; i < ngrid; ++i) {
        double log_f = log(f[i] + F_EPSILON);
        double u, v, mr, mi;
        while ((log_f2 < log_f) && (j < IOP_NTAB - 1)) {
            ++j;
            log_f1  = IOP_tab[j].log_f;
            log_f2  = IOP_tab[j+1].log_f;
            dlog_f  = log_f2 - log_f1;
            mr1     = IOP_tab[j].mr;
            mr2     = IOP_tab[j+1].mr;
            log_mi1 = IOP_tab[j].log_mi;
            log_mi2 = IOP_tab[j+1].log_mi;
        }
        u = (log_f - log_f1) / dlog_f;
        v = 1.0 - u;
        mr = u * mr2 + v * mr1;
        mi = exp(u * log_mi2 + v * log_mi1);
        epsr[i] = mr * mr - mi * mi;
        epsi[i] = 2. * mr * mi;
    }
    return 0;
}   /* H2O_ice_permittivity_Warren_Brandt() */


#ifdef UNIT_TEST

/***********************************************************
* int main(int argc, char **argv)
*
* Purpose:
*   Tests H2O_ice_density() and H2O_ice_permittivity() in
*   this file.  Usage is
*
*   a.out T[C] fmin[GHz] fmax[GHz] num_points
*
*   where num_points is the number of logarithmically-
*   spaced frequency grid points.
*
* Example:
*   $ gcc -D UNIT_TEST h2o_ice.c
*   $ ./a.out -30 0.1 1500 100
************************************************************/


int main(int argc, char **argv)
{
    double T, fmin, fmax, r = 1.0;
    double *f, *epsr, *epsi;
    gridsize_t i;
    int n;

    if (argc < 5) {
        printf("usage: %s T[C] fmin[GHz] fmax[GHz] num_points\n",
            argv[0]);
        return 1;
    }
    T    = atof(argv[1]) + T_STP;
    fmin = atof(argv[2]);
    fmax = atof(argv[3]);
    n    = atoi(argv[4]);
    fmin = fmin < (0.5 * F_EPSILON)  ? (0.5 * F_EPSILON) : fmin;
    if (n > 1)
        r = pow(fmax / fmin, 1. / (double)(n - 1));
    f    = (double*)malloc(n * sizeof(double));
    epsr = (double*)malloc(n * sizeof(double));
    epsi = (double*)malloc(n * sizeof(double));
    printf("# rho(%.8g K) = %.8g cm^-3 (%.8g kg m^-3)\n",
        T,
        H2O_ice_density(T),
        H2O_ice_density(T) * 1.0e6 * MASS_H2O * PCONST_AMU);
    f[0] = fmin;
    for (i = 1; i < n; ++i)
        f[i] = f[i-1] * r;
    H2O_ice_permittivity(epsr, epsi, f, n, T);
    printf("#%13s %14s %14s\n", "f[GHz]", "epsr", "epsi");
    for (i = 0; i < n; ++i)
        printf("%14e %14e %14e\n", f[i], epsr[i], epsi[i]);
    free(f);
    free(epsr);
    free(epsi);
    return 0;
}   /* main() */


/***********************************************************
* int errlog(const int errnum, const int data)
*
* Purpose:
*   dummy errlog() function for unit test
************************************************************/

int errlog(const int errnum, const int data)
{
    printf("errlog(%d, %d)\n", errnum, data);
    return 0;
}   /* errlog() */

/***********************************************************
*int Rayleigh_mol_abscoeff(
*    double *k,
*    const double *f,
*    const gridsize_t ngrid,
*    const double *epsr,
*    const double *epsi,
*    const double rho)
*
* Purpose:
*   dummy Rayleigh_mol_abscoeff() function for unit test
************************************************************/

int Rayleigh_mol_abscoeff(
    double *k,
    const double *f,
    const gridsize_t ngrid,
    const double *epsr,
    const double *epsi,
    const double rho)
{
    return 0;
}   /* Rayleigh_mol_abscoeff() */

#endif /* UNIT_TEST */
