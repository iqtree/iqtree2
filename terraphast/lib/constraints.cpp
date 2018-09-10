#include <terraces/constraints.hpp>

#include "io_utils.hpp"
#include "trees_impl.hpp"
#include "union_find.hpp"

namespace terraces {

std::ostream& operator<<(std::ostream& s, const constraint& c) {
	s << "lca(" << c.left << "," << c.shared << ") < lca(" << c.shared << "," << c.right << ")";
	return s;
}

std::ostream& operator<<(std::ostream& stream, utils::named_output<constraints, name_map> output) {
	auto c = output.entry;
	auto& n = *output.names;
	stream << "lca(" << n[c.left] << "," << n[c.shared] << ") < lca(" << n[c.shared] << ","
	       << n[c.right] << ")";
	return stream;
}

constraints compute_constraints(const std::vector<tree>& trees) {
	constraints result;
	auto num_nodes =
	        (*std::max_element(trees.begin(), trees.end(), [](const tree& a, const tree& b) {
		        return a.size() < b.size();
		})).size();
	std::vector<std::pair<index, index>> outermost_nodes(num_nodes, {none, none});

	for (auto& t : trees) {
		// collect outermost nodes for each subtree (these have lca i)
		foreach_postorder(t, [&](index i) {
			auto node = t[i];
			if (is_leaf(node)) {
				outermost_nodes[i].first = i;
				outermost_nodes[i].second = i;
			} else {
				outermost_nodes[i].first = outermost_nodes[node.lchild()].first;
				outermost_nodes[i].second = outermost_nodes[node.rchild()].second;
			}
		});

		// extract constraints for each edge
		foreach_preorder(t, [&](index i) {
			auto node = t[i];
			if (!is_leaf(node)) {
				auto lchild = node.lchild();
				auto rchild = node.rchild();
				auto lnode = t[lchild];
				auto rnode = t[rchild];
				// taxon of leftmost descendant of i
				auto i1 = t[outermost_nodes[i].first].taxon();
				// taxon of rightmost descendant of lchild of i
				auto i2 = t[outermost_nodes[lchild].second].taxon();
				// taxon of leftmost descendant of rchild of i
				auto i3 = t[outermost_nodes[rchild].first].taxon();
				// taxon of rightmost descendant of i
				auto i4 = t[outermost_nodes[i].second].taxon();

				// if the left edge is an inner edge
				if (!is_leaf(lnode)) {
					result.emplace_back(i2, i1, i4);
				}
				// if the right edge is an inner edge
				if (!is_leaf(rnode)) {
					result.emplace_back(i3, i4, i1);
				}
			}
		});
	}

	return result;
}

index deduplicate_constraints(constraints& in_c) {
	for (auto& c : in_c) {
		c = {std::min(c.left, c.shared), std::max(c.left, c.shared), c.right};
	}
	std::sort(in_c.begin(), in_c.end(), [](constraint a, constraint b) {
		return std::tie(a.left, a.shared, a.right) < std::tie(b.left, b.shared, b.right);
	});
	auto it = std::unique(in_c.begin(), in_c.end());
	index count = (index)std::distance(it, in_c.end());
	in_c.erase(it, in_c.end());
	return count;
}

} // namespace terraces
