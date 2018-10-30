#include <terraces/bitmatrix.hpp>

#include <cassert>

namespace terraces {

bitmatrix::bitmatrix(index rows, index cols) : m_rows{rows}, m_cols{cols}, m_vec(rows * cols) {}

index bitmatrix::rows() const { return m_rows; }
index bitmatrix::cols() const { return m_cols; }

bool bitmatrix::get(index row, index col) const {
	assert(row < m_rows && col < m_cols);
	return m_vec[row * m_cols + col];
}

void bitmatrix::set(index row, index col, bool val) {
	assert(row < m_rows && col < m_cols);
	m_vec[row * m_cols + col] = val;
}

void bitmatrix::row_or(index in1, index in2, index out) {
	for (index i = 0; i < cols(); ++i) {
		set(out, i, get(in1, i) || get(in2, i));
	}
}

bitmatrix bitmatrix::get_cols(const std::vector<std::size_t>& cols) const {
	assert(cols.size() <= this->cols());
	auto ret = bitmatrix{rows(), cols.size()};
	for (auto i = std::size_t{}; i < rows(); ++i) {
		for (auto j = std::size_t{}; j < cols.size(); ++j) {
			ret.set(i, j, get(i, cols[j]));
		}
	}
	return ret;
}

bool bitmatrix::operator==(const bitmatrix& other) const {
	return other.rows() == rows() && other.cols() == cols() && other.m_vec == m_vec;
}

bool bitmatrix::operator!=(const bitmatrix& other) const { return !(*this == other); }

} // namespace terraces
