/*

BOOSTER: BOOtstrap Support by TransfER: 
BOOSTER is an alternative method to compute bootstrap branch supports 
in large trees. It uses transfer distance between bipartitions, instead
of perfect match.

Copyright (C) 2017 Frederic Lemoine, Jean-Baka Domelevo Entfellner, Olivier Gascuel

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef _PARSIMONY_UTILS
#define _PARSIMONY_UTILS

#include "tree.h"
#include "hashmap.h"
#include "sort.h"
#include "stats.h"

Tree* gen_random_tree(Tree *tree);
/**
   This function precomputes the esperence of the expected number of parsimony steps
   implied by a bipartition under the hypothesis that the tree is random.
   In Input:
   - The max depth
   - The number of taxa
   - A pointer to a 2D array (given by precompute_steps_probability(int max_depth, int nb_tax)):
      * First dimension : depth
      * Second dimension : steps
      * value : probability of the step at a given depth
   In output : An array with :
   - the depth in index
   - the expected Number of random parsimony steps  
*/

#endif
