/*
 *  sine_fast.h
 *  pbcrlibm
 *
 *  Created by Catherine Daramy on Thu Mar 04 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "crlibm.h"
#include "crlibm_private.h"

/* some useful vlues */
const double twotopihi = 0x3FE45F306DC9C883;
const double twotopilo = 0xBC86B01EC5418056;

const double piotwohi = 0x3FF921FB54442D18;
const double piotwolo = 0x3C91A62633145C06;

db_number const poly_sin_fast[9] =
{
{{0x3FF00000,0x00000000}},
{{0xBFC55555,0x55555555}},
{{0x3F811111,0x11111111}},
{{0xBF2A01A0,0x1A01A01A}},
{{0x3EC71DE3,0xA556C734}},
{{0xBE5AE645,0x67F544E4}},
{{0x3DE61246,0x13A86D09}},
{{0xBD6AE7F3,0xE733B81F}},
{{0x3CE952C7,0x7030AD4A}}
};

