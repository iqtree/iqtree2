#include <terraces/parser.hpp>

#include <algorithm>
#include <array>
#include <stack>
#include <stdexcept>
#include <utility>
#include <vector>

#include <terraces/errors.hpp>

#include "trees_impl.hpp"
#include "utils.hpp"

namespace terraces {

// non-public parsing-utilities:
namespace parsing {

enum class token_type { lparen, rparen, name, seperator, eof };

struct token {
	token_type type;
	std::string name;

	token(token_type type, std::string name = "") : type{type}, name{std::move(name)} {}
};

struct parser_state {
	index parent;
	index self;

	parser_state(index parent, index self) : parent{parent}, self{self} {}
};

using parser_stack = std::stack<parser_state, std::vector<parser_state>>;

template <typename Iterator>
token next_token(Iterator& it, Iterator end) {
	constexpr static std::array<char, 3> special_tokens{{'(', ')', ','}};
	it = utils::skip_ws(it, end);
	if (it == end) {
		return {token_type::eof};
	}
	if (*it != '\'') {
		switch (*it) {
		case '(':
			++it;
			return {token_type::lparen};
		case ')':
			++it;
			return {token_type::rparen};
		case ',':
			++it;
			return {token_type::seperator};
		}
	} else {
		++it;
		const auto name_begin = it;
		it = std::find(it, end, '\'');
		const auto name_end = it;
		utils::ensure<bad_input_error>(name_end != end,
		                               bad_input_error_type::nwk_mismatched_quotes,
		                               std::string{name_begin, name_end});
		++it;
		return {token_type::name, {name_begin, name_end}};
	}
	const auto name_begin = it;
	it = std::find_first_of(it, end, special_tokens.begin(), special_tokens.end());
	const auto name_end = utils::reverse_skip_ws(name_begin, std::find(name_begin, it, ':'));
	return {token_type::name, {name_begin, name_end}};
}

template <typename NameCallback>
tree parse_nwk_impl(const std::string& input, NameCallback cb) {
	auto ret = tree{};

	auto stack = parsing::parser_stack{};

	auto it = input.begin();
	const auto end = input.end();

	bool unrooted = false;

	ret.emplace_back(none, none, none, none);
	auto state = parsing::parser_state{none, 0};

	for (auto token = parsing::next_token(it, end); token.type != parsing::token_type::eof;
	     token = parsing::next_token(it, end)) {
		switch (token.type) {
		case parsing::token_type::lparen: {
			const auto parent = state.self;
			const auto self = ret.size();
			// Inner node names must come AFTER the closing parentheses!
			utils::ensure<bad_input_error>(ret[state.self].taxon() == none,
			                               bad_input_error_type::nwk_malformed);
			stack.push(state);
			state = parsing::parser_state{parent, self};
			ret.emplace_back(parent, none, none, none);
			ret[parent].lchild() = self;
			break;
		}
		case parsing::token_type::seperator: {
			const auto parent = state.parent;
			if (ret[parent].rchild() != none) {
				// (*,old,new)root, state = {root,old}
				// -> (*,(old,new)aux)root
				utils::ensure<bad_input_error>(
				        parent == 0, bad_input_error_type::nwk_multifurcating);
				unrooted = true;
				auto old_node = state.self;
				auto aux_node = ret.size();
				auto new_node = ret.size() + 1;
				assert(ret[0].rchild() == state.self);
				ret.emplace_back(0, old_node, new_node, none);
				ret.emplace_back(aux_node, none, none, none);
				ret[0].rchild() = aux_node;
				ret[old_node].parent() = aux_node;
				stack.push({0, aux_node});
				state.parent = aux_node;
				state.self = new_node;
			} else {
				state.self = ret.size();
				ret.emplace_back(parent, none, none, none);
				ret[parent].rchild() = state.self;
			}
			// no need to update state as the tree is binary to
			// begin with, which means that we will now go up a level
			break;
		}
		case parsing::token_type::rparen: {
			utils::ensure<bad_input_error>(
			        not stack.empty(),
			        bad_input_error_type::nwk_mismatched_parentheses);
			utils::ensure<bad_input_error>(ret[state.parent].rchild() != none,
			                               bad_input_error_type::nwk_malformed);
			state = stack.top();
			stack.pop();
			break;
		}
		case parsing::token_type::name: {
			cb(ret[state.self], token.name);
			break;
		}
		case parsing::token_type::eof:
		default: { throw std::logic_error{"dafuq?"}; }
		}
	}
	if (unrooted) {
		utils::ensure<bad_input_error>(!stack.empty(),
		                               bad_input_error_type::nwk_mismatched_parentheses);
		stack.pop();
	}
	utils::ensure<bad_input_error>(stack.empty(),
	                               bad_input_error_type::nwk_mismatched_parentheses);
	return ret;
}

} // namespace parsing

tree parse_nwk(const std::string& input, const index_map& taxa) {
	std::vector<bool> found_taxon(taxa.size(), false);
	return parsing::parse_nwk_impl(input, [&](node& n, const std::string& name) {
		if (is_leaf(n)) {
			auto it = taxa.find(name);
			utils::ensure<bad_input_error>(
			        it != taxa.end(), bad_input_error_type::nwk_taxon_unknown, name);
			auto taxon_id = (*it).second;
			utils::ensure<bad_input_error>(!found_taxon[taxon_id],
			                               bad_input_error_type::nwk_taxon_duplicate,
			                               name);
			found_taxon[taxon_id] = true;
			n.taxon() = taxon_id;
		}
	});
}

named_tree parse_new_nwk(const std::string& input) {
	name_map names;
	index_map indices;
	auto t = parsing::parse_nwk_impl(input, [&](node& n, const std::string& name) {
		if (is_leaf(n)) {
			auto ret = indices.insert({name, names.size()});
			utils::ensure<bad_input_error>(
			        ret.second, bad_input_error_type::nwk_taxon_duplicate, name);
			n.taxon() = names.size();
			names.emplace_back(name);
		}
	});
	return {t, names, indices};
}

occurrence_data parse_bitmatrix(std::istream& input) {
	index cols{};
	index rows{};
	input >> rows >> cols >> std::ws;

	bitmatrix mat{rows, cols};
	name_map names;
	index_map indices;

	std::string line{};
	while (std::getline(input, line)) {
		if (line.empty()) {
			continue;
		}
		auto it = line.begin();
		auto end = line.end();
		auto taxon_id = names.size();

		// fill matrix
		for (index i = 0; i < cols; ++i) {
			it = utils::skip_ws(it, end);
			utils::ensure<bad_input_error>(
			        it != end, bad_input_error_type::bitmatrix_size_invalid, line);
			auto c = *it++;
			utils::ensure<bad_input_error>(c == '1' || c == '0',
			                               bad_input_error_type::bitmatrix_malformed);
			if (c == '1') {
				mat.set(taxon_id, i, true);
			}
		}

		// read taxon name
		it = utils::skip_ws(it, end);
		utils::ensure<bad_input_error>(it != end,
		                               bad_input_error_type::bitmatrix_name_empty, line);
		auto taxon_name = std::string{it, end};
		auto was_inserted = indices.insert({taxon_name, names.size()}).second;
		utils::ensure<bad_input_error>(
		        was_inserted, bad_input_error_type::bitmatrix_name_duplicate, taxon_name);
		names.emplace_back(std::move(taxon_name));
	}
	utils::ensure<bad_input_error>(rows == names.size(),
	                               bad_input_error_type::bitmatrix_size_invalid);

	return {mat, names, indices};
}

} // namespace terraces
