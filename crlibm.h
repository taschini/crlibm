/*
 * Author  : Defour David, Catherine Daramy
 * Contact : David.Defour@ens-lyon.fr, catherine_daramy@ens-lyon.fr
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
 */

#ifndef CRLIBM_H
#define CRLIBM_H


#ifdef CRLIBM_TYPECPU_X86
#include <fpu_control.h>
#ifndef __setfpucw
#define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
#endif 
#endif

/* An init function which sets FPU flags when needed (mostly on Intel
   architectures with default double extended) */
extern void crlibm_init(void);

/* Compute the cosine with correct rounding */
extern double cos_rn(double); /* to nearest  */
extern double cos_rd(double); /* toward -inf */ 
extern double cos_ru(double); /* toward +inf */ 

/* Compute the sine with correct rounding */
extern double sin_rn(double); /* to nearest  */
extern double sin_rd(double); /* toward -inf */ 
extern double sin_ru(double); /* toward +inf */ 

/* Compute the tangent with correct rounding */
extern double tan_rn(double); /* to nearest  */
extern double tan_rd(double); /* toward -inf */ 
extern double tan_ru(double); /* toward +inf */ 

/* Compute the tangent with correct rounding */
extern double cotan_rn(double); /* to nearest  */
extern double cotan_rd(double); /* toward -inf */ 
extern double cotan_ru(double); /* toward +inf */ 

#if 0 /* temporarily disabled */
/* Compute the arctangent with correct rounding */
extern double atan_rn(double); /* to nearest  */
extern double atan_rd(double); /* toward -inf */ 
extern double atan_ru(double); /* toward +inf */ 
#endif

/* Compute the exponential with correct rounding */
extern double exp_rn(double); /* to nearest  */
extern double exp_rd(double); /* toward -inf */ 
extern double exp_ru(double); /* toward +inf */ 
#define exp_rz exp_rd         /* toward zero */ 

extern double exp2_rn(double); /* to nearest  */
extern double exp2_rd(double); /* toward -inf */ 
extern double exp2_ru(double); /* toward +inf */ 

/* Compute the logarithm with correct rounding */
extern double log_rn(double); /* to nearest  */
extern double log_rd(double); /* toward -inf */ 
extern double log_ru(double); /* toward +inf */ 

extern double log2_rn(double); /* to nearest  */
extern double log2_rd(double); /* toward -inf */ 
extern double log2_ru(double); /* toward +inf */ 

extern double log10_rn(double); /* to nearest  */
extern double log10_rd(double); /* toward -inf */ 
extern double log10_ru(double); /* toward +inf */ 



#endif /* ifdef CRLIBM_H*/
