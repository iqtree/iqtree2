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

#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `nrows' to `rows' and `ncols' to `cols'. In addition, memory is allocated for `data' (each element of 
|	the matrix `data' is a NxsDiscreteDatum object, which can do its own initialization).
*/
NxsDiscreteMatrix::NxsDiscreteMatrix(
  unsigned rows,	/* number of taxa */
  unsigned cols)	/* number of characters */
	{
	nrows = rows;
	ncols = cols;

	data = new NxsDiscreteDatum*[nrows];
	for (unsigned i = 0; i < nrows; i++)
		data[i] = new NxsDiscreteDatum[ncols];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes memory allocated in the constructor for data member `data'.
*/
NxsDiscreteMatrix::~NxsDiscreteMatrix()
	{
	if (data != NULL)
		{
		for (unsigned i = 0; i < nrows; i++)
			delete [] data[i];
		delete [] data;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for `nAddRows' additional rows and updates the variable nrows. Data already stored in `data' is not
|	destroyed; the newly-allocated rows are added at the bottom of the existing matrix.
*/
void NxsDiscreteMatrix::AddRows(
  unsigned nAddRows)	/* the number of additional rows to allocate */
	{
	unsigned new_nrows = nrows + nAddRows;

	// Allocate a matrix big enough to hold all of the existing data
	// as well as the new rows.
	//
	NxsDiscreteDatum **new_data = new NxsDiscreteDatum*[new_nrows];

	// Copy existing data to the new matrix.
	//
	unsigned i;
	for (i = 0; i < nrows; i++)
		new_data[i] = data[i];

	// Let data now point to the new data matrix
	//
	delete [] data;
	data = new_data;

	// Create data elements for the newly added rows
	//
	for (i = nrows; i < new_nrows; i++)
		data[i] = new NxsDiscreteDatum[ncols];

	nrows = new_nrows;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds state directly to the NxsDiscreteDatum object at data[i][j]. Assumes `data' is non-NULL, `i' is in the range 
|	[0..nrows), and `j' is in the range [0..ncols). The `value' argument is assumed to be either zero or a positive 
|	integer. Calls private member function AddState to do the real work; look at the documentation for that function 
|	for additional details.
*/
void NxsDiscreteMatrix::AddState(
  unsigned i,		/* the (0-offset) index of the taxon in question */
  unsigned j,		/* the (0-offset) index of the character in question */
  unsigned value)	/* the state to be added */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);
	assert(value >= 0);

	AddState(data[i][j], value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds an additional state to the array states of `d'. If `states' is NULL, allocates memory for two integers and 
|	assigns 1 to the first and `value' to the second. If `states' is non-NULL, allocates a new int array long enough to 
|	hold states already present plus the new one being added here, then deletes the old `states' array. Assumes that we 
|	are not trying to set either the missing state or the gap state here; the functions SetMissing or SetGap, 
|	respectively, should be used for those purposes. Also assumes that we do not want to overwrite the state. This 
|	function always adds states to those already present; use SetState to overwrite the state.
*/
void NxsDiscreteMatrix::AddState(
  NxsDiscreteDatum &d,	/* the NxsDiscreteDatum object affected */
  unsigned value)		/* the additional state to be added */
	{
	unsigned oldns = GetNumStates(d);
	unsigned k, newlen;

	unsigned *tmp = d.states;

	if (IsMissing(d))
		{
		d.states = new unsigned[2];
		d.states[0] = 1;
		d.states[1] = value;
		}

	else if (IsGap(d))
		{
		d.states = new unsigned[2];
		d.states[0] = 1;
		d.states[1] = value;
		}

	else if (oldns == 1)
		{
		d.states = new unsigned[4];
		d.states[0] = 2;
		d.states[1] = tmp[1];
		d.states[2] = value;
		d.states[3] = 0;	// Assume not polymorphic unless told otherwise.
		}

	else
		{
		newlen = oldns + 3;
		d.states = new unsigned[newlen];
		d.states[0] = oldns + 1;
		for (k = 1; k < oldns + 1; k++)
		d.states[k] = tmp[k];
		d.states[newlen - 2] = value;
		d.states[newlen - 1] = 0;	// Assume not polymorphic unless told otherwise.
		}

	if (tmp != NULL)
		delete [] tmp;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets state of taxon `i' and character `j' to state of first taxon for character `j'. Assumes `i' is in the range 
|	[0..nrows) and `j' is in the range [0..ncols). Also assumes `data' is non-NULL. Calls private function CopyFrom
|	to do the actual work.
*/
void NxsDiscreteMatrix::CopyStatesFromFirstTaxon(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	data[i][j].CopyFrom(data[0][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a dump of the current contents of the data matrix stored in the variable `data'. Translates missing data
|	elements to the '?' character and gap states to '-', otherwise, calls GetState to provide the representation.
*/
void NxsDiscreteMatrix::DebugSaveMatrix(
  ostream &out,			/* the stream on which to dump the matrix contents */
  unsigned colwidth)	/* the width of a data matrix column in characters */
	{
	out << endl;
	out << "nrows = " << nrows << endl;
	out << "ncols = " << ncols << endl;
	for (unsigned i = 0; i < nrows; i++)
		{
		for (unsigned j = 0; j < ncols; j++)
			{
			if (IsMissing(i, j))
				out << setw(colwidth) << '?';
			else if (IsGap(i, j))
				out << setw(colwidth) << '-';
			else
				out << setw(colwidth) << GetState(i, j);
			}
		out << endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Duplicates columns `startCol' to `endCol' in row `row' of the matrix. If additional storage is needed to accommodate
|	the duplication, this is done automatically through the use of the AddRows method. Note that `count' includes the 
|	row already present, so if `count' is 10, then 9 more rows will actually be added to the matrix to make a total of 
|	10 identical rows. The parameters `startCol' and `endCol' default to 0 and `ncols', so if duplication of the entire
|	row is needed, these need not be explicitly specified in the call to DuplicateRow. Return value is number of 
|	additional rows allocated to matrix (0 if no rows needed to be allocated). Assumes `data' is non-NULL, `row' is in
|	the range [[0..`nrows'), `startCol' is in the range [0..`ncols'), and `endCol' is either UINT_MAX, in which case it
|	is reset to `ncols' - 1, or is in the range (`startCol'..`ncols').
*/
unsigned NxsDiscreteMatrix::DuplicateRow(
  unsigned row,			/* the row to be duplicated */
  unsigned count,		/* the total number of copies needed */
  unsigned startCol,	/* the starting column (inclusive) in the range of columns to be duplicated */
  unsigned endCol)		/* the ending column (inclusive) in the range of columns to be duplicated */
	{
	assert(data != NULL);
	assert(row >= 0);
	assert(row < nrows);
	assert(startCol >= 0);
	assert(startCol < ncols);
	if (endCol == UINT_MAX)
		endCol = ncols - 1;
	assert(endCol > startCol);
	assert(endCol < ncols);

	// Expand matrix (if necessary) to accommodate additional rows.
	//
	unsigned nNewRows = 0;
	if (row + count > nrows)
		{
		nNewRows = row + count - nrows;
		AddRows(nNewRows);
		}

	for (unsigned i = 1; i < count; i++)
		{
		for (unsigned col = startCol; col <= endCol; col++)
			data[row+i][col] = data[row][col];
		}

	return nNewRows;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all cells of `data', setting `data' to NULL, and resets `nrows' and `ncols' to 0.
*/
void NxsDiscreteMatrix::Flush()
	{
	if (data != NULL)
		{
		for (unsigned i = 0; i < nrows; i++)
			delete [] data[i];
		delete [] data;
		}

	nrows	= 0;
	ncols	= 0;
	data	= NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes that `data' is non-NULL, `i' is in the range [0..`nrows') and `j' is in the range [0..`ncols'). Returns 
|	reference to the NxsDiscreteDatum object at row `i', column `j' of matrix.
*/
NxsDiscreteDatum &NxsDiscreteMatrix::GetDiscreteDatum(
  unsigned i,	/* the row of the matrix */
  unsigned j)	/* the column of the matrix */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	return data[i][j];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of states for taxon `i' and character `j'. Assumes `data' is non-NULL, `i' is in the range 
|	[0..`nrows'), and `j' is in the range [0..`ncols'). Calls private member function GetNumStates to do the actual
|	work.
*/
unsigned NxsDiscreteMatrix::GetNumStates(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	return GetNumStates(data[i][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns total number of states assigned to `d'. Returns 0 for both gap and missing states.
*/
unsigned NxsDiscreteMatrix::GetNumStates(
  NxsDiscreteDatum &d)	/* the datum in question */
	{
	if (d.states == NULL)
		return 0;

	return d.states[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of states for character `j' over all taxa. Note: this function is rather slow, as it must walk 
|	through each taxon for the specified character, adding the states encountered to a set, then finally returning the
|	size of the set. Thus, if this function is called often, it would be advisable to initialize an array using this 
|	function, then refer to the array subsequently. Assumes `j' is in the range [0..`ncols') and `data' is non-NULL.
|	Includes all taxa (i.e. there is no mechanism here for treating some taxa as deleted for a particular analysis).
|	Missing and gap states are ignored.
*/
unsigned NxsDiscreteMatrix::GetObsNumStates(
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	// Create a set object to hold all states seen for all taxa for character j
	//
	set< unsigned, less<unsigned> > stateset;

	for (unsigned i = 0; i < nrows; i++)
		{
		NxsDiscreteDatum &d = data[i][j];
		unsigned ns = GetNumStates(d);
		if (ns == 0)
			continue;
		for (unsigned k = 0; k < ns; k++)
			stateset.insert(GetState(d, k));
		}

	return stateset.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `k'th state possessed by taxon `i' and character `j'. This taxon-character combination will have more 
|	than one state if there is ambiguity or polymorphism. Assumes that `i' is in the range [0..`nrows') and `j' is in 
|	the range [0..`ncols'). Also assumes that at least one state is present (i.e., not the gap or missing state). Use 
|	the function GetNumStates to determine the  number of states present. Assumes `k' is in the range [0..ns), where ns
|	is the value returned by GetNumStates.
*/
unsigned NxsDiscreteMatrix::GetState(
  unsigned i,	/* the row of the matrix */
  unsigned j,	/* the column of the matrix */
  unsigned k)	/* the state to return */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	return GetState(data[i][j], k);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the internal unsigned representation of the state stored in `d' at position `k' of the array `d.states'. 
|	Assumes that the state is not the missing or gap state. Use IsMissing and IsGap prior to calling this function to 
|	ensure this function will succeed. Assumes that `k' is in the range [ 0 .. `d.states'[0]).
*/
unsigned NxsDiscreteMatrix::GetState(
  NxsDiscreteDatum &d,	/* the datum in question */
  unsigned k)			/* the number of the state */
	{
	assert(!IsMissing(d));
	assert(!IsGap(d));
	assert(k >= 0);
	assert(k < d.states[0]);

	return d.states[k + 1];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns 1 if the state for taxon `i', character `j', is set to the gap symbol, 0 otherwise. Assumes `data' is 
|	non-NULL, `i' is in the range [0..`nrows') and `j' is in the range [0..`ncols').
*/
bool NxsDiscreteMatrix::IsGap(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	return IsGap(data[i][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the gap state is stored, otherwise returns false. Note: returns false if this datum represents 
|	missing data (often the gap state is equated with missing data, but the distinction is made here).
*/
bool NxsDiscreteMatrix::IsGap(
  NxsDiscreteDatum &d)	/* the datum in question */
	{
	if (d.states == NULL || d.states[0] > 0)
		return 0;
	else
		return 1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns 1 if the state for taxon `i', character `j', is set to the missing data symbol, 0 otherwise. Assumes `i' is 
|	in the range [0..`nrows') and `j' is in the range [0..`ncols').
*/
bool NxsDiscreteMatrix::IsMissing(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	return IsMissing(data[i][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the missing state is stored, false otherwise. Note that this function returns false if the gap state
|	is stored (often the gap state is equated with missing data, but the distinction is maintained here).
*/
bool NxsDiscreteMatrix::IsMissing(
  NxsDiscreteDatum &d)	/* the datum in question */
	{
	if (d.states == NULL)
		return 1;
	else
		return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns 1 if character `j' is polymorphic in taxon `i', 0 otherwise. Assumes `data' is non-NULL, `i' is in the 
|	range [0..`nrows') and `j' is in the range [0..`ncols').
*/
bool NxsDiscreteMatrix::IsPolymorphic(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	return IsPolymorphic(data[i][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the number of states is greater than 1 and polymorphism has been specified. Returns false if the 
|	state stored is the missing state, the gap state, or if the number of states is 1.
*/
bool NxsDiscreteMatrix::IsPolymorphic(
  NxsDiscreteDatum &d)	/* the datum in question */
	{
	if (d.states == NULL || d.states[0] < 2)
		return 0;

	int nstates = d.states[0];
	int ncells = nstates + 2;
	return (bool)(d.states[ncells - 1] > 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all cells of `data' and reallocates memory to create a new matrix object with `nrows' = `rows' and `ncols' 
|	= `cols'. Assumes `rows' and `cols' are both greater than 0.
*/
void NxsDiscreteMatrix::Reset(
  unsigned rows,	/* the new number of rows (taxa) */
  unsigned cols)	/* the new number of columns (characters) */
	{
	unsigned i;
	assert(rows > 0);
	assert(cols > 0);

	// Delete what is there now
	//
	if (data != NULL)
		{
		for (i = 0; i < nrows; i++)
			delete [] data[i];
		delete [] data;
		}

	nrows = rows;
	ncols = cols;

	// Create new data matrix
	//
	data = new NxsDiscreteDatum*[nrows];
	for (i = 0; i < nrows; i++)
		data[i] = new NxsDiscreteDatum[ncols];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets state stored at `data[i][j]' to the gap state. Assumes `i' is in the range [0..`nrows') and `j' is in the 
|	range [0..`ncols'). Calls the private SetGap member function to do the actual work.
*/
void NxsDiscreteMatrix::SetGap(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	SetGap(data[i][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns the gap state to `d', erasing any previously stored information. The gap state is designated internally as 
|	a states array one element long, with the single element set to the value 0.
*/
void NxsDiscreteMatrix::SetGap(
  NxsDiscreteDatum &d)	/* the datum in question */
	{
	if (d.states != NULL)
		delete [] d.states;
	d.states = new unsigned[1];
	d.states[0] = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets state stored at `data[i][j]' to the missing state. Assumes `data' is non-NULL, `i' is in the range [0..`nrows')
|	and `j' is in the range [0..`ncols'). Calls the private member function SetMissing to do the actual work.
*/
void NxsDiscreteMatrix::SetMissing(
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	SetMissing(data[i][j]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns the missing state to `d', erasing any previously stored information. The missing state is stored internally
|	as a NULL value for the states array.
*/
void NxsDiscreteMatrix::SetMissing(
  NxsDiscreteDatum &d)	/* the datum in question */
	{
	if (d.states != NULL)
		delete [] d.states;
	d.states = NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specify 1 for `value' if taxon at row `i' is polymorphic at character in column `j', 0 for `value' if uncertain 
|	which state applies. Sets polymorphism state of taxon `i' and character `j' to `value'. Assumes `data' is non-NULL,
|	`i' is in the range [0..`nrows') and `j' is in the range [0..`ncols'). Also assumes that the number of states 
|	stored is greater than 1. Calls private member function SetPolymorphic to do the actual work.
*/
void NxsDiscreteMatrix::SetPolymorphic(
  unsigned i,		/* the (0-offset) index of the taxon in question */
  unsigned j,		/* the (0-offset) index of the character in question */
  unsigned value)	/* specify either 0 or 1, where 0 means ambiguity and 1 means polymorphism */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);
	assert(value == 0 || value == 1);

	SetPolymorphic(data[i][j], value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the polymorphism cell (last cell in `d.states') to `value'. Warning: has no effect if there are fewer than 2 
|	states stored!
*/
void NxsDiscreteMatrix::SetPolymorphic(
  NxsDiscreteDatum &d,	/* the datum in question */
  unsigned value)		/* specify 1 if polymorphic, 0 if uncertain */
	{
	if (d.states == NULL || d.states[0] < 2)
		return;

	int nstates = d.states[0];
	int ncells = nstates + 2;
	d.states[ncells - 1] = value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets state of taxon `i' and character `j' to `value'. Assumes `data' is non-NULL, `i' is in the range [0..`nrows') 
|	and `j' is in the range [0..`ncols'). Assumes that this function will not be called if there is missing data or the 
|	state is the gap state, in which case the functions SetMissing or SetGap, respectively, should be called instead.
|	Calls the private member function SetState to do the actual work.
*/
void NxsDiscreteMatrix::SetState(
  unsigned i,		/* the (0-offset) index of the taxon in question */
  unsigned j,		/* the (0-offset) index of the character in question */
  unsigned value)	/* the value to assign for this state */
	{
	assert(i >= 0);
	assert(i < nrows);
	assert(j >= 0);
	assert(j < ncols);
	assert(data != NULL);

	SetState(data[i][j], value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns `value' to the 2nd cell in `d.states' (1st cell in `d.states' array is set to 1 to indicate that there is 
|	only one state). Warning: if already one or more states (including the gap state) are assigned to `d', they will 
|	be forgotten. Use the function AddState if you want to preserve states already stored in `d'. Assumes state being 
|	set is not the missing state nor the gap state; use SetMissing or SetGap, respectively, for that.
*/
void NxsDiscreteMatrix::SetState(
  NxsDiscreteDatum &d,	/* the datum in question */
  unsigned value)		/* the value to assign for the state */
	{
	if (d.states != NULL)
		delete [] d.states;
	d.states = new unsigned[2];
	d.states[0] = 1;
	d.states[1] = value;
	}

