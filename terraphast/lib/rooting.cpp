#include "trees_impl.hpp"
#include "utils.hpp"
#include "validation.hpp"
#include <sstream>
#include <stack>
#include <terraces/rooting.hpp>

namespace terraces {

namespace {

bool is_rchild(const tree& t, index node) { return t[t[node].parent()].rchild() == node; }

void copy_subtree(const tree& src, tree& dst, index sub_root) {
	std::stack<index> stack;
	stack.push(sub_root);
	while (!stack.empty()) {
		auto idx = stack.top();
		auto node = src[idx];
		stack.pop();
		if (is_leaf(node)) {
			dst[idx].lchild() = none;
			dst[idx].rchild() = none;
			dst[idx].taxon() = node.taxon();
		} else {
			dst[idx].lchild() = node.lchild();
			dst[idx].rchild() = node.rchild();
			dst[idx].taxon() = none;
			dst[node.lchild()].parent() = idx;
			dst[node.rchild()].parent() = idx;
			stack.push(node.lchild());
			stack.push(node.rchild());
		}
	}
}

void copy_reversed(const tree& src, tree& dst, index child, index cur, index parent) {
	auto coming_from_right = is_rchild(src, child);
	auto opposite_idx = src[cur].child(!coming_from_right);
	dst[cur].child(coming_from_right) = parent;
	dst[cur].child(!coming_from_right) = opposite_idx;
	dst[parent].parent() = cur;
	dst[opposite_idx].parent() = cur;
	copy_subtree(src, dst, opposite_idx);
}

void copy_reversed_final(const tree& src, tree& dst, index child, index cur, index root_rchild) {
	auto coming_from_right = is_rchild(src, child);
	auto opposite_idx = src[cur].child(!coming_from_right);
	dst[cur].child(coming_from_right) = root_rchild;
	dst[cur].child(!coming_from_right) = opposite_idx;
	dst[root_rchild].parent() = cur;
	dst[opposite_idx].parent() = cur;
	copy_subtree(src, dst, opposite_idx);
	copy_subtree(src, dst, root_rchild);
}

void insert_root(tree& t, index lchild, index rchild, bool reverse) {
	t[0].parent() = none;
	t[0].taxon() = none;
	t[0].child(reverse) = lchild;
	t[0].child(!reverse) = rchild;
	t[lchild].parent() = 0;
	t[rchild].parent() = 0;
}

} // anonymous namespace

std::vector<bool> root_split(const tree& t, index num_leaves) {
	std::vector<bool> split(num_leaves);
	auto root = t[0];
	foreach_preorder(t,
	                 [&](index i) {
		                 auto node = t[i];
		                 if (is_leaf(node)) {
			                 split[node.taxon()] = true;
		                 }
		         },
	                 root.rchild());
	return split;
}

tree reroot_at_node(const tree& t, index node_idx) {
	utils::ensure<std::invalid_argument>(node_idx != 0, "can't reroot at the root");
	terraces::check_rooted_tree(t);
	// if the node is already below the root, we don't reroot.
	if (t[node_idx].parent() == 0) {
		return t;
	}

	tree out(t.size());

	// everything below node_idx is unchanged
	copy_subtree(t, out, node_idx);

	// the root is placed between node_idx and its parent
	// (keeping the orientation of node_idx wrt its parent)
	bool is_node_rchild = is_rchild(t, node_idx);
	insert_root(out, node_idx, t[node_idx].parent(), is_node_rchild);

	auto prev = node_idx;
	auto cur = t[node_idx].parent();
	auto next = t[cur].parent();
	while (next != 0) {
		copy_reversed(t, out, prev, cur, next);
		prev = utils::exchange(cur, utils::exchange(next, t[next].parent()));
	}

	// everything on the opposite side of the original root is unchanged
	auto root_rchild = is_rchild(t, cur);
	auto opposite_idx = t[0].child(!root_rchild);
	copy_reversed_final(t, out, prev, cur, opposite_idx);

	check_rooted_tree(out);

	assert(is_isomorphic_unrooted(t, out));

	return out;
}

void reroot_at_taxon_inplace(tree& t, index comp_taxon) {
	// identify node corresponding to comp_taxon
	index root_leaf = none;
	for (index i = 0; i < t.size(); ++i) {
		if (t[i].taxon() == comp_taxon) {
			assert(root_leaf == none);
			root_leaf = i;
		}
	}
	assert(root_leaf != none && "The tree doesn't contain the given taxon");
	terraces::check_rooted_tree(t);

	// first: swap at inner nodes s.t. root_leaf is the rightmost leaf
	index cur = root_leaf;
	index p = t[cur].parent();
	while (cur != 0) {
		// if cur lies in the left subtree, swap it to the right
		if (t[p].lchild() == cur) {
			std::swap(t[p].lchild(), t[p].rchild());
		}
		cur = utils::exchange(p, t[p].parent());
	}
	// second: move the root down right until we meet root_leaf
	auto& li = t[0].lchild();
	auto& ri = t[0].rchild();
	while (ri != root_leaf) {
		auto& ln = t[li];
		auto& rn = t[ri];
		auto& rli = rn.lchild();
		auto& rri = rn.rchild();
		auto& rrn = t[rri];
		std::swap(rrn.parent(), ln.parent());
		std::tie(li, ri, rli, rri) = std::make_tuple(ri, rri, li, rli);
	}
}

} // namespace terraces
