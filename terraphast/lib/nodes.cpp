#include <terraces/trees.hpp>

#include <ostream>

namespace terraces {

std::ostream& operator<<(std::ostream& s, const node& n) {
	auto print_index = [&s](index i) {
		if (i == none) {
			s << "none";
		} else {
			s << i;
		}
	};
	s << '(';
	print_index(n.data[0]);
	s << ", ";
	print_index(n.data[1]);
	s << ", ";
	print_index(n.data[2]);
	s << ')';
	return s;
}

} // namespace terraces
