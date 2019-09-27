#include "../lib/trees_impl.hpp"
#include <algorithm>
#include <fstream>
#include <terraces/advanced.hpp>
#include <terraces/bitmatrix.hpp>
#include <terraces/errors.hpp>
#include <terraces/parser.hpp>
#include <terraces/terraces.h>

missingData* initializeMissingData(size_t numberOfSpecies, size_t numberOfPartitions,
                                   const char** speciesNames) {
	auto data = new missingData;
	data->numberOfSpecies = numberOfSpecies;
	data->numberOfPartitions = numberOfPartitions;
	data->allocatedNameArray = false; // TODO What is this entry?
	data->speciesNames = speciesNames;
	data->missingDataMatrix = new unsigned char[numberOfSpecies * numberOfPartitions]();
	return data;
}

void freeMissingData(missingData* m) {
	delete[] m->missingDataMatrix;
	delete m;
}

void setDataMatrix(missingData* m, size_t speciesNumber, size_t partitionNumber,
                   unsigned char value) {
	m->missingDataMatrix[speciesNumber * m->numberOfPartitions + partitionNumber] = value;
}

unsigned char getDataMatrix(const missingData* m, size_t speciesNumber, size_t partitionNumber) {
	return m->missingDataMatrix[speciesNumber * m->numberOfPartitions + partitionNumber];
}

void copyDataMatrix(const unsigned char* matrix, missingData* m) {
	std::copy_n(matrix, m->numberOfSpecies * m->numberOfPartitions, m->missingDataMatrix);
}

CHECK_RESULT int terraceAnalysis(missingData* m, const char* newickTreeString, const int ta_outspec,
                                 const char* allTreesOnTerraceFile, mpz_t terraceSize) {
	// check ta_outspec
	auto detect = bool(ta_outspec & TA_DETECT);
	auto count = bool(ta_outspec & TA_COUNT);
	auto enumerate = bool(ta_outspec & TA_ENUMERATE);
	auto compress = bool(ta_outspec & TA_ENUMERATE_COMPRESS);
	auto force_comprehensive = bool(ta_outspec & TA_UPPER_BOUND);
	bool invalid1 = detect && (count || enumerate); // cannot detect and count at the same time
	bool invalid2 = compress && !enumerate;         // cannot compress if we don't enumerate
	if (invalid1 || invalid2) {
		return TERRACE_FLAG_CONFLICT_ERROR;
	}

	// check input sizes
	if (m->numberOfPartitions < 2) {
		return TERRACE_NUM_PARTITIONS_ERROR;
	}
	if (m->numberOfSpecies < 4) {
		return TERRACE_NUM_SPECIES_ERROR;
	}

	// copy missing data matrix
	terraces::bitmatrix matrix{m->numberOfSpecies, m->numberOfPartitions};
	for (size_t row = 0; row < m->numberOfSpecies; ++row) {
		size_t rowcount = 0;
		for (size_t col = 0; col < m->numberOfPartitions; ++col) {
			auto val = m->missingDataMatrix[row * m->numberOfPartitions + col];
			if (val != 0 && val != 1) {
				return TERRACE_MATRIX_ERROR;
			}
			matrix.set(row, col, val);
			rowcount += val;
		}
		if (rowcount == 0) {
			return TERRACE_SPECIES_WITHOUT_PARTITION_ERROR;
		}
	}

	// copy names
	terraces::name_map names;
	terraces::index_map name_index;
	for (size_t spec_i = 0; spec_i < m->numberOfSpecies; ++spec_i) {
		names.emplace_back(m->speciesNames[spec_i]);
		if (!name_index.insert({names.back(), spec_i}).second) {
			return TERRACE_SPECIES_ERROR;
		}
	}

	// parse newick tree
	terraces::tree tree;
	try {
		tree = terraces::parse_nwk(newickTreeString, name_index);
	} catch (const terraces::bad_input_error& err) {
		switch (err.type()) {
		case terraces::bad_input_error_type::nwk_multifurcating:
			return TERRACE_TREE_NOT_BINARY_ERROR;
		case terraces::bad_input_error_type::nwk_taxon_duplicate:
			return TERRACE_SPECIES_ERROR;
		default:
			return TERRACE_NEWICK_ERROR;
		}
	}
	if (terraces::num_leaves_from_nodes(tree.size()) != m->numberOfSpecies) {
		return TERRACE_SPECIES_ERROR;
	}

	// prepare data
	if (force_comprehensive) {
		matrix = terraces::maximum_comprehensive_columnset(matrix);
	}

	terraces::supertree_data data;
	try {
		data = terraces::create_supertree_data(tree, matrix);
	} catch (const terraces::bad_input_error&) {
		return TERRACE_INTERNAL_ERROR;
	} catch (const terraces::no_usable_root_error&) {
		return TERRACE_NO_ROOT_SPECIES_ERROR;
	}

	// enumerate terrace
	if (detect) {
		auto lb = terraces::fast_count_terrace(data);
		mpz_set_ui(terraceSize, lb);
	} else if (count && !enumerate) {
		try {
			auto size = terraces::count_terrace_bigint(data);
			mpz_set(terraceSize, size.value().get_mpz_t());
		} catch (const terraces::tree_count_overflow_error&) {
			return TERRACE_SPLIT_COUNT_OVERFLOW_ERROR;
		}
	} else {
		auto ofs = std::ofstream{allTreesOnTerraceFile};
		if (not ofs.is_open()) {
			return TERRACE_OUTPUT_FILE_ERROR;
		}
		mpz_class size;
		try {
			if (compress) {
				size = terraces::print_terrace_compressed(data, names, ofs).value();
			} else {
				size = terraces::print_terrace(data, names, ofs).value();
			}
			if (count) {
				mpz_set(terraceSize, size.get_mpz_t());
			}
		} catch (std::ifstream::failure&) {
			return TERRACE_OUTPUT_FILE_ERROR;
		} catch (const terraces::tree_count_overflow_error&) {
			return TERRACE_SPLIT_COUNT_OVERFLOW_ERROR;
		}
	}
	return TERRACE_SUCCESS;
}
