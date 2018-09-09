#ifndef SUPERTREE_VARIANTS_DEBUG_HPP
#define SUPERTREE_VARIANTS_DEBUG_HPP

#include <algorithm>
#include <ostream>
#include <vector>

#include "constraints_impl.hpp"
#include "io_utils.hpp"
#include "supertree_variants.hpp"
#include "union_find_debug.hpp"

namespace terraces {
namespace debug {
namespace variants {

template <typename Result>
struct stack_state {
	index current_bip;
	index max_bip;
	bool right;
	Result intermediate;

	stack_state(Result intermediate)
	        : current_bip{0}, max_bip{0}, right{false}, intermediate{intermediate} {}
};

template <typename Callback>
class stack_state_decorator : public Callback {
public:
	using result_type = typename Callback::result_type;

private:
	std::vector<stack_state<result_type>> m_stack;

public:
	template <typename... Args>
	stack_state_decorator(Args&&... args) : Callback{std::forward<Args>(args)...} {}

	result_type begin_iteration(const bipartitions& bip_it, const bitvector& c_occ,
	                            const constraints& constraints) {
		auto result = Callback::begin_iteration(bip_it, c_occ, constraints);
		m_stack.emplace_back(result);
		m_stack.back().current_bip = bip_it.begin_bip();
		m_stack.back().max_bip = bip_it.end_bip();
		return result;
	}

	void step_iteration(const bipartitions& bip_it, index bip) {
		Callback::step_iteration(bip_it, bip);
		m_stack.back().current_bip = bip;
	}

	void finish_iteration() {
		m_stack.pop_back();
		Callback::finish_iteration();
	}

	void right_subcall() {
		Callback::right_subcall();
		m_stack.back().right = true;
	}

	void left_subcall() {
		Callback::left_subcall();
		m_stack.back().right = false;
	}

	result_type accumulate(result_type acc, result_type val) {
		auto result = Callback::accumulate(acc, val);
		m_stack.back().intermediate = result;
		return result;
	}

	const std::vector<stack_state<result_type>>& stack() const { return m_stack; }
};

template <typename Callback>
class logging_decorator : public Callback {
public:
	using result_type = typename Callback::result_type;

private:
	std::ostream& m_output;
	std::size_t m_depth;
	bool m_first_iteration;
	const name_map& m_names;

	std::ostream& output() { return m_output << std::string(m_depth, '\t'); }

public:
	logging_decorator(Callback&& cb, std::ostream& output, const name_map& names)
	        : Callback{std::move(cb)}, m_output{output}, m_depth{0}, m_first_iteration{},
	          m_names{names} {}

	void enter(const ranked_bitvector& leaves) {
		Callback::enter(leaves);
		output() << "<itr>\n";
		++m_depth;
		output() << "<leafs val=\"{" << utils::as_comma_separated_output(leaves, m_names)
		         << "}\" />\n";
	}

	result_type begin_iteration(const bipartitions& bip_it, const bitvector& c_occ,
	                            const constraints& c) {
		output() << "<constraints val=\"{"
		         << utils::as_comma_separated_output(c_occ, c, m_names) << "}\" />\n";
		output() << "<sets val=\"[";
		bool set_first = true;
		auto& leaves = bip_it.leaves();
		for (auto rep : debug::union_find_iterable_sets{bip_it.sets()}) {
			if (not set_first) {
				m_output << ",";
			}
			m_output << "{";
			bool el_first = true;
			for (auto el : debug::union_find_iterable_set{bip_it.sets(), rep}) {
				if (not el_first) {
					m_output << ",";
				}
				m_output << m_names[leaves.select(el)];
				el_first = false;
			}
			m_output << "}";
			set_first = false;
		}
		m_output << "]\" />\n";
		m_first_iteration = true;

		return Callback::begin_iteration(bip_it, c_occ, c);
	}

	void step_iteration(const bipartitions& bip_it, index bip) {
		Callback::step_iteration(bip_it, bip);
		if (not m_first_iteration) {
			--m_depth;
			output() << "</bipartition>\n";
		}
		auto subleaf_sets = bip_it.get_both_sets_unsafe(bip);
		output() << "<bipartition l=\"{"
		         << utils::as_comma_separated_output(subleaf_sets.first, m_names)
		         << "}\" r=\"{";
		m_output << utils::as_comma_separated_output(subleaf_sets.second, m_names)
		         << "}\" />\n";
		++m_depth;
		m_first_iteration = false;
	}

	void finish_iteration() {
		--m_depth;
		output() << "</bipartition>\n";
		Callback::finish_iteration();
	}

	result_type exit(result_type result) {
		output() << "<result val=\"" << result << "\" />\n";
		--m_depth;
		output() << "</itr>\n";
		return Callback::exit(result);
	}
};

} // namespace variants
} // namespace debug
} // namespace terraces

#endif // SUPERTREE_VARIANTS_DEBUG_HPP
