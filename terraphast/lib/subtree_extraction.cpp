#include "trees_impl.hpp"
#include "utils.hpp"
#include <terraces/errors.hpp>
#include <terraces/subtree_extraction.hpp>

using std::vector;
using std::stack;
using std::pair;
using std::make_pair;

namespace terraces {

std::pair<bitmatrix, std::vector<index>> compute_node_occ(const tree& t, const bitmatrix& occ) {
	auto num_nodes = num_nodes_from_leaves(occ.rows());
	auto num_sites = occ.cols();
	utils::ensure<bad_input_error>(t.size() == num_nodes,
	                               bad_input_error_type::tree_mismatching_size);
	check_rooted_tree(t);
	auto node_occ = bitmatrix{t.size(), occ.cols()};
	std::vector<index> num_leaves_per_site(occ.cols(), 0);

	// compute occurrences on inner nodes: bitwise or of the children
	foreach_postorder(t, [&](index i) {
		auto node = t[i];
		if (is_leaf(node)) {
			// copy data from taxon occurrence matrix
			utils::ensure<bad_input_error>(node.taxon() != none,
			                               bad_input_error_type::tree_unnamed_leaf);
			for (index site = 0; site < num_sites; ++site) {
				auto has_leaf = occ.get(node.taxon(), site);
				node_occ.set(i, site, has_leaf);
				num_leaves_per_site[site] += has_leaf;
			}
		} else {
			node_occ.row_or(node.lchild(), node.rchild(), i);
		}
	});
	return {std::move(node_occ), std::move(num_leaves_per_site)};
}

index induced_lca(const tree& t, const bitmatrix& node_occ, index column) {
	index lca = 0;
	while (!is_leaf(t[lca])) {
		auto present = [&](index node) { return node_occ.get(node, column); };
		assert(present(lca));
		auto node = t[lca];
		if (present(node.lchild()) && present(node.rchild())) {
			return lca;
		} else {
			lca = present(node.lchild()) ? node.lchild() : node.rchild();
		}
	}
	return lca;
}

tree subtree(const tree& t, const bitmatrix& node_occ,
             const std::vector<index>& num_leaves_per_site, index site) {
	auto root = induced_lca(t, node_occ, site);
	if (is_leaf(t[root])) {
		// tree containing only a single leaf
		return {{none, none, none, t[root].taxon()}};
	}

	auto present = [&](index node) { return node_occ.get(node, site); };
	assert(present(t[root].lchild()) && present(t[root].rchild()));

	tree out_tree;
	out_tree.reserve(num_nodes_from_leaves(num_leaves_per_site[site]));
	out_tree.emplace_back(); // root node

	stack<index> boundary;
	auto callback = [&](index i) {
		auto node = t[i];
		bool leaf_occ = is_leaf(node) && present(i);
		bool inner_occ = !is_leaf(node) && present(node.lchild()) && present(node.rchild());

		if (leaf_occ || (inner_occ && i != root)) {
			// fires if the tree is trivial (i.e. only one edge!)
			// this can only happen with sites for which only one
			// species has data, which we should have caught earlier.
			assert(!boundary.empty());
			auto parent = boundary.top();
			out_tree.emplace_back(parent, none, none, node.taxon());
			if (out_tree[parent].lchild() == none) {
				out_tree[parent].lchild() = out_tree.size() - 1;
			} else {
				assert(out_tree[parent].rchild() == none);
				out_tree[parent].rchild() = out_tree.size() - 1;
				boundary.pop();
			}
		}
		if (inner_occ) {
			boundary.push(out_tree.size() - 1);
		}
	};
	foreach_preorder(t, callback, root);

	return out_tree;
}

std::vector<tree> subtrees(const tree& t, const bitmatrix& occ) {
	auto num_sites = occ.cols();
	// TODO naming of compute_node_occ
	const auto node_occ_pair = compute_node_occ(t, occ);
	const auto& node_occ = node_occ_pair.first;
	const auto& num_leaves_per_site = node_occ_pair.second;

	// collect leaves and inner nodes: bitwise and of the children
	vector<tree> out_trees;

	for (index site = 0; site < num_sites; ++site) {
		out_trees.push_back(subtree(t, node_occ, num_leaves_per_site, site));
	}

	return out_trees;
}

} // namespace terraces
