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

#include "phylotree.h"
#include "tinatree.h"
#include "parsmultistate.h"
#include "alignment.h"

void doParsMultiState(Params &params) {
	cout << "Here\n";
    Alignment alignment(params.aln_file, params.sequence_type, params.intype);
    TinaTree tree;
    tree.readTree(params.user_file, params.is_rooted);
	tree.setAlignment(&alignment);
	tree.drawTree(cout);
	cout << "Parsimony score is: " << tree.computeParsimonyScore() << endl;
	cout << "Parsimony score ver2 is: " << tree.computeParsimony() << endl;
	//tree.printParsimonyStates();
}
