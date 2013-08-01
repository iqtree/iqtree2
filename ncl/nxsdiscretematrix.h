//	Copyright (C) 1999-2003 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library) version 2.0.
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef NCL_NXSDISCRETEMATRIX_H
#define NCL_NXSDISCRETEMATRIX_H

#include <limits.h>

/*----------------------------------------------------------------------------------------------------------------------
|	Class providing storage for the discrete data types (dna, rna, nucleotide, standard, and protein) inside a DATA or 
|	CHARACTERS block. This class is also used to store the data for an ALLELES block. Maintains a matrix in which each 
|	cell is an object of the class NxsDiscreteDatum. NxsDiscreteDatum stores the state for a particular combination of 
|	taxon and character as an integer. Ordinarily, there will be a single state recorded for each taxon-character 
|	combination, but exceptions exist if there is polymorphism for a taxon-character combination, or if there is 
|	uncertainty about the state (e.g., in dna data, the data file might have contained an R or Y entry). Please consult 
|	the documentation for the NxsDiscreteDatum class for the details about how states are stored. For data stored in an 
|	ALLELES block, rows of the matrix correspond to individuals and columns to loci. Each NxsDiscreteDatum must 
|	therefore store information about both genes at a single locus for a single individual in the case of diploid data.
|	To do this, two macros HIWORD and LOWORD are used to divide up the unsigned value into two words. A maximum of 255 
|	distinct allelic forms can be accommodated by this scheme, assuming at minimum a 32-bit architecture. Because it is
|	not known in advance how many rows are going to be necessary, The NxsDiscreteMatrix class provides the AddRows 
|	method, which expands the number of rows allocated for the matrix while preserving data already stored. 
*/
class NxsDiscreteMatrix
	{
	friend class NxsCharactersBlock;
	friend class NxsAllelesBlock;

	public:

							NxsDiscreteMatrix(unsigned rows, unsigned cols);
		virtual				~NxsDiscreteMatrix();

		void				AddRows(unsigned nAddRows);
		void				AddState(unsigned i, unsigned j, unsigned value);
		void				CopyStatesFromFirstTaxon(unsigned i, unsigned j);
		void				DebugSaveMatrix(ostream &out, unsigned colwidth = 12);
		unsigned			DuplicateRow(unsigned row, unsigned count, unsigned startCol = 0, unsigned endCol = UINT_MAX);
		void				Flush();
		unsigned			GetState(unsigned i, unsigned j, unsigned k = 0);
		unsigned			GetNumStates(unsigned i, unsigned j);
		unsigned			GetObsNumStates(unsigned j);
		bool				IsGap(unsigned i, unsigned j);
		bool				IsMissing(unsigned i, unsigned j);
		bool				IsPolymorphic(unsigned i, unsigned j);
		void				Reset(unsigned rows, unsigned cols);
		void				SetGap(unsigned i, unsigned j);
		void				SetMissing(unsigned i, unsigned j);
		void				SetPolymorphic(unsigned i, unsigned j, unsigned value = 1);
		void				SetState(unsigned i, unsigned j, unsigned value);

	private:

		unsigned			nrows;	/* number of rows (taxa) in the data matrix */
		unsigned			ncols;	/* number of columns (characters) in the data matrix */
		NxsDiscreteDatum	**data;	/* storage for the data */

		void				AddState(NxsDiscreteDatum &d, unsigned value);
		bool				IsGap(NxsDiscreteDatum &d);
		bool				IsMissing(NxsDiscreteDatum &d);
		bool				IsPolymorphic(NxsDiscreteDatum &d);
		NxsDiscreteDatum	&GetDiscreteDatum(unsigned i, unsigned j);
		unsigned			GetNumStates(NxsDiscreteDatum &d);
		unsigned			GetState(NxsDiscreteDatum &d, unsigned k = 0);
		void				SetGap(NxsDiscreteDatum &d);
		void				SetMissing(NxsDiscreteDatum &d);
		void				SetPolymorphic(NxsDiscreteDatum &d, unsigned value);
		void				SetState(NxsDiscreteDatum &d, unsigned value);
	};

typedef NxsDiscreteMatrix DiscreteMatrix;


#endif
