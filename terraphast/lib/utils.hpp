
#ifndef TERRACES_UTILS_HPP
#define TERRACES_UTILS_HPP

#include <cctype>
#include <iterator>
#include <utility>

namespace terraces {
namespace utils {

template <typename Iterator>
Iterator skip_ws(Iterator it, Iterator last) {
	while (it != last and std::isspace(*it)) {
		++it;
	}
	return it;
}

template <typename Iterator>
Iterator reverse_skip_ws(Iterator first, Iterator last) {
	while (first != last) {
		if (std::isspace(*std::prev(last))) {
			--last;
		} else {
			break;
		}
	}
	return last;
}

template <typename Exception, typename... Args>
void ensure(bool b, Args&&... args) {
	if (not b) {
		throw Exception{std::forward<Args>(args)...};
	}
}

} // namespace utils
} // namespace terraces

#endif // TERRACES_UTILS_HPP
