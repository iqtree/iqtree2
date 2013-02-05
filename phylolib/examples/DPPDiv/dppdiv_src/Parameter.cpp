/* 
 * DPPDiv version 1.0b source code (git: 9c0ac3d2258f89827cfe9ba2b5038f0f656b82c1)
 * Copyright 2009-2011
 * Tracy Heath(1,2,3) (NSF postdoctoral fellowship in biological informatics DBI-0805631)
 * Mark Holder(1)
 * John Huelsenbeck(2)
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 * DPPDiv is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck
 */

#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_exchangeability.h"
#include "Parameter_expcalib.h"
#include "Parameter_rate.h"
#include "Parameter_shape.h"
#include "Parameter_tree.h"
#include "Parameter_cphyperp.h"
#include "Parameter_treescale.h"
#include "Parameter_speciaton.h"

using namespace std;

Parameter::Parameter(MbRandom *rp, Model *mp) {

	ranPtr = rp;
	modelPtr = mp;
}

Parameter& Parameter::operator=(Parameter &p) {

	if (this != &p){
		ranPtr = p.ranPtr;
		name   = p.name;
		
		{
			Basefreq *thatDerivedPtr = dynamic_cast<Basefreq *>(&p);
			Basefreq *thisDerivedPtr = dynamic_cast<Basefreq *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}

		{
			Exchangeability *thatDerivedPtr = dynamic_cast<Exchangeability *>(&p);
			Exchangeability *thisDerivedPtr = dynamic_cast<Exchangeability *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}

		{
			Shape *thatDerivedPtr = dynamic_cast<Shape *>(&p);
			Shape *thisDerivedPtr = dynamic_cast<Shape *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}
		
		{
			Tree *thatDerivedPtr = dynamic_cast<Tree *>(&p);
			Tree *thisDerivedPtr = dynamic_cast<Tree *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}

		{
			NodeRate *thatDerivedPtr = dynamic_cast<NodeRate *>(&p);
			NodeRate *thisDerivedPtr = dynamic_cast<NodeRate *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}
			
		{
			Cphyperp *thatDerivedPtr = dynamic_cast<Cphyperp *>(&p);
			Cphyperp *thisDerivedPtr = dynamic_cast<Cphyperp *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}
		
		{
			Treescale *thatDerivedPtr = dynamic_cast<Treescale *>(&p);
			Treescale *thisDerivedPtr = dynamic_cast<Treescale *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}
		
		{
			Speciation *thatDerivedPtr = dynamic_cast<Speciation *>(&p);
			Speciation *thisDerivedPtr = dynamic_cast<Speciation *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}
		
		{
			ExpCalib *thatDerivedPtr = dynamic_cast<ExpCalib *>(&p);
			ExpCalib *thisDerivedPtr = dynamic_cast<ExpCalib *>(this);
			if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 ){
				thisDerivedPtr->clone( *thatDerivedPtr );
				goto exitOperator;
			}
		}
		exitOperator:
			;
	}
		
	return *this;
}