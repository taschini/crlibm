/*
 *  sinh.h
 *  
 *
 *  Created by Catherine Daramy on Fri Mar 14 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */
#include "crlibm.h"
#include "crlibm_private.h"

static const scs bigconst =
/* ~7.104759e+02 */ 
{{0x000002c6, 0x1e747dcf, 0x30fb6207, 0x36f7e921, 
0x381a56ef, 0x3bf0039c, 0x0622ac09, 0x0f5a8642},
DB_ONE,   0,   1 } 
,
mediumconst =
/* ~5.500000e+01 */ 
{{0x00000037, 0x00000000, 0x00000000, 0x00000000, 
0x00000000, 0x00000000, 0x00000000, 0x00000000},
DB_ONE,   0,   1 } 
,
miniconst =
/* ~9.094947e-13 */ 
{{0x00100000, 0x00000000, 0x00000000, 0x00000000, 
0x00000000, 0x00000000, 0x00000000, 0x00000000},
DB_ONE,  -2,   1 };


#define miniconst_ptr (scs_ptr)&miniconst
#define mediumconst_ptr (scs_ptr)&mediumconst
#define bigconst_ptr (scs_ptr)&bigconst