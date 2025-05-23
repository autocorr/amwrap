/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* oneline.c                   S. Paine rev. 2018 February 16
*
* Line catalog and partition sum table for oneline, a test
* molecule with a one-line spectrum, used for program
* validation.
************************************************************/

#include "oneline.h"
#include "molecules.h"

/*
 * Oneline has only one isotopologue, with a very small mass to
 * give the same half width at half maximum for Doppler and
 * pressure-broadened lines at 500 GHz, 500 mbar, and 296 K.
 */

const double oneline_abundance_tab[] = {
    1.0,
    ABUN_ONELINE,
};

const double oneline_mass_tab[] = {
    MASS_ONELINE,
    MASS_ONELINE,
};

const double oneline_Tref = 296.0;

/*
 * The line catalog.  See am_types.h for the definition of the
 * cat_entry_t structure.
 */

const cat_entry_t oneline_cat[] = {
{    500.000000, 1.000000e-020,    0.000000, 1.0000e-003f, 2.0000e-003f, 0.70f,  0.0000e+000f,  1 },
};

const int oneline_num_lines = (int)(sizeof(oneline_cat) / sizeof(cat_entry_t));

/*
 * Partition sum table.  In each row, column 0 is the temperature T [K].
 * Oneline has only one isotopologue, with Q in column 1
 */

const double oneline_Qtab[] = {
    200.0,  1.886937,
    201.0,  1.887467,
    202.0,  1.887991,
    203.0,  1.888511,
    204.0,  1.889026,
    205.0,  1.889536,
    206.0,  1.890042,
    207.0,  1.890543,
    208.0,  1.891040,
    209.0,  1.891532,
    210.0,  1.892019,
    211.0,  1.892502,
    212.0,  1.892981,
    213.0,  1.893456,
    214.0,  1.893926,
    215.0,  1.894393,
    216.0,  1.894855,
    217.0,  1.895313,
    218.0,  1.895767,
    219.0,  1.896218,
    220.0,  1.896664,
    221.0,  1.897107,
    222.0,  1.897546,
    223.0,  1.897981,
    224.0,  1.898412,
    225.0,  1.898840,
    226.0,  1.899265,
    227.0,  1.899685,
    228.0,  1.900103,
    229.0,  1.900516,
    230.0,  1.900927,
    231.0,  1.901334,
    232.0,  1.901737,
    233.0,  1.902138,
    234.0,  1.902535,
    235.0,  1.902929,
    236.0,  1.903320,
    237.0,  1.903707,
    238.0,  1.904092,
    239.0,  1.904473,
    240.0,  1.904852,
    241.0,  1.905227,
    242.0,  1.905600,
    243.0,  1.905969,
    244.0,  1.906336,
    245.0,  1.906700,
    246.0,  1.907061,
    247.0,  1.907419,
    248.0,  1.907775,
    249.0,  1.908128,
    250.0,  1.908478,
    251.0,  1.908825,
    252.0,  1.909170,
    253.0,  1.909512,
    254.0,  1.909852,
    255.0,  1.910189,
    256.0,  1.910524,
    257.0,  1.910856,
    258.0,  1.911186,
    259.0,  1.911513,
    260.0,  1.911838,
    261.0,  1.912160,
    262.0,  1.912480,
    263.0,  1.912798,
    264.0,  1.913114,
    265.0,  1.913427,
    266.0,  1.913738,
    267.0,  1.914047,
    268.0,  1.914353,
    269.0,  1.914658,
    270.0,  1.914960,
    271.0,  1.915260,
    272.0,  1.915558,
    273.0,  1.915854,
    274.0,  1.916148,
    275.0,  1.916440,
    276.0,  1.916729,
    277.0,  1.917017,
    278.0,  1.917303,
    279.0,  1.917587,
    280.0,  1.917869,
    281.0,  1.918149,
    282.0,  1.918427,
    283.0,  1.918703,
    284.0,  1.918977,
    285.0,  1.919250,
    286.0,  1.919521,
    287.0,  1.919789,
    288.0,  1.920056,
    289.0,  1.920322,
    290.0,  1.920585,
    291.0,  1.920847,
    292.0,  1.921107,
    293.0,  1.921366,
    294.0,  1.921622,
    295.0,  1.921877,
    296.0,  1.922131,
    297.0,  1.922382,
    298.0,  1.922633,
    299.0,  1.922881,
    300.0,  1.923128,
    301.0,  1.923373,
    302.0,  1.923617,
    303.0,  1.923859,
    304.0,  1.924100,
    305.0,  1.924339,
    306.0,  1.924577,
    307.0,  1.924813,
    308.0,  1.925048,
    309.0,  1.925281,
    310.0,  1.925513,
    311.0,  1.925743,
    312.0,  1.925972,
    313.0,  1.926200,
    314.0,  1.926426,
    315.0,  1.926651,
    316.0,  1.926874,
    317.0,  1.927096,
    318.0,  1.927317,
    319.0,  1.927536,
    320.0,  1.927754,
    321.0,  1.927971,
    322.0,  1.928187,
    323.0,  1.928401,
    324.0,  1.928614,
    325.0,  1.928825,
    326.0,  1.929036,
    327.0,  1.929245,
    328.0,  1.929453,
    329.0,  1.929660,
    330.0,  1.929865,
    331.0,  1.930069,
    332.0,  1.930272,
    333.0,  1.930474,
    334.0,  1.930675,
    335.0,  1.930875,
    336.0,  1.931073,
    337.0,  1.931271,
    338.0,  1.931467,
    339.0,  1.931662,
    340.0,  1.931856,
    341.0,  1.932049,
    342.0,  1.932241,
    343.0,  1.932431,
    344.0,  1.932621,
    345.0,  1.932810,
    346.0,  1.932997,
    347.0,  1.933184,
    348.0,  1.933369,
    349.0,  1.933553,
    350.0,  1.933737,
    351.0,  1.933919,
    352.0,  1.934101,
    353.0,  1.934281,
    354.0,  1.934460,
    355.0,  1.934639,
    356.0,  1.934816,
    357.0,  1.934993,
    358.0,  1.935168,
    359.0,  1.935343,
    360.0,  1.935517,
    361.0,  1.935690,
    362.0,  1.935861,
    363.0,  1.936032,
    364.0,  1.936202,
    365.0,  1.936371,
    366.0,  1.936540,
    367.0,  1.936707,
    368.0,  1.936873,
    369.0,  1.937039,
    370.0,  1.937204,
    371.0,  1.937368,
    372.0,  1.937531,
    373.0,  1.937693,
    374.0,  1.937854,
    375.0,  1.938014,
    376.0,  1.938174,
    377.0,  1.938333,
    378.0,  1.938491,
    379.0,  1.938648,
    380.0,  1.938805,
    381.0,  1.938960,
    382.0,  1.939115,
    383.0,  1.939269,
    384.0,  1.939422,
    385.0,  1.939575,
    386.0,  1.939727,
    387.0,  1.939877,
    388.0,  1.940028,
    389.0,  1.940177,
    390.0,  1.940326,
    391.0,  1.940474,
    392.0,  1.940621,
    393.0,  1.940768,
    394.0,  1.940913,
    395.0,  1.941059,
    396.0,  1.941203,
    397.0,  1.941347,
    398.0,  1.941490,
    399.0,  1.941632,
    400.0,  1.941773,
};

const int oneline_Qtab_cols = 2;
const int oneline_Qtab_rows =
    (int)(sizeof(oneline_Qtab) / (2 * sizeof(double)));
