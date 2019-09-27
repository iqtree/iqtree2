/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "tree/phylotree.h"
#include "tree/tinatree.h"
#include "parsmultistate.h"
#include "alignment/alignment.h"
#include "tree/parstree.h"

void doParsMultiState(Params &params) {
    Alignment alignment(params.aln_file, params.sequence_type, params.intype, "");
    alignment.orderPatternByNumChars(PAT_VARIANT);
    ParsTree pars_tree;
    pars_tree.readTree(params.user_file, params.is_rooted);
    if (pars_tree.rooted)
        pars_tree.convertToUnrooted();
    pars_tree.setAlignment(&alignment);
    pars_tree.initCostMatrix(CM_LINEAR);
    pars_tree.setParsimonyKernel(params.SSE);
    pars_tree.initializeAllPartialPars();
    int total_length = pars_tree.computeParsimony();
    cout << "total length: " << total_length << endl;
    pars_tree.initCostMatrix(CM_UNIFORM);
    int pars_score = pars_tree.computeParsimony();
    cout.unsetf(ios::fixed);
    cout.precision(6);
    cout << "mean length: " << double(total_length)/pars_score << endl;
    cout << "Parsimony score is: " << pars_score << endl;
	//cout << "Parsimony score ver2 is: " << tree.computeParsimony() << endl;
	//tree.printParsimonyStates();
}
