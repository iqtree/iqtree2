#ifndef TERRACES_TREES_HPP
#define TERRACES_TREES_HPP

#include <array>
#include <cassert>
#include <cstdint>
#include <iosfwd>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "definitions.hpp"

namespace terraces {

/**
 * Represents a node of a rooted tree.
 * Its three neighbor indices are stored in \ref data and
 * in the rooted case can be accessed using \ref lchild, \ref rchild
 * and \ref parent.
 */
struct node {
	node() : data{{none, none, none, none}} {}
	node(index parent, index left, index right, index taxon)
	        : data{{parent, left, right, taxon}} {}
	/* data[0]: parent
	     * data[1]: left child
	     * data[2]: right child
	 * data[3]: taxon index */
	std::array<index, 4> data = {{none, none, none, none}};

	index parent() const { return data[0]; }
	index& parent() { return data[0]; }

	index lchild() const { return data[1]; }
	index& lchild() { return data[1]; }

	index rchild() const { return data[2]; }
	index& rchild() { return data[2]; }

	index child(bool right) const { return data[1u + bool(right)]; }
	index& child(bool right) { return data[1u + bool(right)]; }

	index taxon() const { return data[3]; }
	index& taxon() { return data[3]; }

	bool operator==(const node& o) const { return data == o.data; }

	bool operator!=(const node& o) const { return !(o == *this); }
};

/** A tree is represented by its nodes. The root of a tree is stored at index 0. */
using tree = std::vector<node>;

/** Stores the name for every node (or leaf) of a tree. */
using name_map = std::vector<std::string>;

/** Helper struct for Newick tree output. */
struct newick_t {
	const tree* t;
	const name_map* names;
};

/** Helper function for Newick tree output. */
inline newick_t as_newick(const tree& t, const name_map& names) { return {&t, &names}; }

/**
 * Prints a Newick tree to the given output stream.
 * @see as_newick
 */
std::ostream& operator<<(std::ostream& s, newick_t tree_pair);

/** Maps the name of a species to its index in the tree. */
using index_map = std::unordered_map<std::string, index>;

} // namespace terraces

#endif // TERRACES_TREES_HPP
