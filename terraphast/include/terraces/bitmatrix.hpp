#ifndef TERRACES_BITMATRIX_HPP
#define TERRACES_BITMATRIX_HPP

#include <vector>

#include "trees.hpp"

namespace terraces {

/** A (memory-wise) compact bitmatrix. */
class bitmatrix {
public:
	/** Constructs a bitmatrix with \p rows rows and \p cols columns. */
	bitmatrix(index rows, index cols);

	/** @returns the number of rows. */
	index rows() const;
	/** @returns the number of columns. */
	index cols() const;

	/** Returns the entry at cordinates (row, col). */
	bool get(index row, index col) const;
	/** Sets the entry at coordinates (row, col). */
	void set(index row, index col, bool val);
	/** Writes the bit rows \p in1 | \p in2 to row \p out. */
	void row_or(index in1, index in2, index out);

	/** Returns a bitmatrix containing only the given columns. */
	bitmatrix get_cols(const std::vector<std::size_t>& cols) const;

	/** Tests for equality with another bitmatrix. */
	bool operator==(const bitmatrix& other) const;
	/** Tests for inequality with another bitmatrix. */
	bool operator!=(const bitmatrix& other) const;

private:
	index m_rows;
	index m_cols;
	std::vector<bool> m_vec; // yes, it's evil
};
}

#endif // TERRACES_BITMATRIX_HPP
