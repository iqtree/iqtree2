#include <iostream>
#include <random>
#include <string>

int main(int argc, char** argv) {
	if (argc != 4) {
		std::cerr << "Usage: " << argv[0]
		          << " <number-of-species> <number-of_sites> <probabilty-of-1>\n";
		return 1;
	}
	const auto num_species = std::stoul(argv[1]);
	const auto num_sites = std::stoul(argv[2]);
	const auto prob = std::stod(argv[3]);

	auto gen = std::default_random_engine{std::random_device{}()};
	auto dist_p = std::uniform_real_distribution<double>{0.0, 1.0};

	const auto root_species =
	        std::uniform_int_distribution<unsigned long>{0u, num_species}(gen);

	std::cout << num_species << ' ' << num_sites << '\n';
	for (auto i = 0ul; i < num_species; ++i) {
		if (i == root_species) {
			for (auto j = 0ul; j < num_sites; ++j) {
				std::cout << "1 ";
			}
		} else {
			for (auto j = 0ul; j < num_sites; ++j) {
				std::cout << (dist_p(gen) < prob) << ' ';
			}
		}
		std::cout << 's' << i << '\n';
	}
}
