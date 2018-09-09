#include <catch.hpp>

#include <iostream>

#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>

#include "../lib/validation.hpp"

namespace terraces {
namespace tests {

TEST_CASE("is_isomorphic_unrooted (simple true)", "[validation]") {
	auto fst = parse_new_nwk("(1,(2,3));");
	reroot_at_taxon_inplace(fst.tree, fst.indices.at("1"));
	auto snd = parse_nwk("(2,(1,3));", fst.indices);
	reroot_at_taxon_inplace(snd, fst.indices.at("3"));
	CHECK(is_isomorphic_unrooted(fst.tree, snd));
}

TEST_CASE("is_isomorphic_unrooted (simple false)", "[validation]") {
	auto fst = parse_new_nwk("(1,(2,(3,(4,5))));");
	reroot_at_taxon_inplace(fst.tree, fst.indices.at("1"));
	auto snd = parse_nwk("(2,((1,4),(3,5)));", fst.indices);
	reroot_at_taxon_inplace(snd, fst.indices.at("3"));
	CHECK(!is_isomorphic_unrooted(fst.tree, snd));
}

TEST_CASE("is_isomorphic_unrooted (complex)", "[validation]") {
	auto fst = parse_new_nwk("((((s2,s4),((s13,s1),s7)),s3),s5);");
	reroot_at_taxon_inplace(fst.tree, fst.indices.at("s2"));
	auto snd = parse_nwk("((s13,((s2,s7),(s4,(s5,s3)))),s1);", fst.indices);
	reroot_at_taxon_inplace(snd, fst.indices.at("s7"));
	auto trd = parse_nwk("((s13,(s7,((s4,s2),(s5,s3)))),s1);", fst.indices);
	reroot_at_taxon_inplace(trd, fst.indices.at("s13"));
	auto fth = parse_nwk("((s13,(((s3,s5),s4),(s2,s7))),s1);", fst.indices);
	reroot_at_taxon_inplace(fth, fst.indices.at("s5"));
	CHECK(!is_isomorphic_unrooted(fst.tree, snd));
	CHECK(is_isomorphic_unrooted(fst.tree, trd));
	CHECK(!is_isomorphic_unrooted(snd, trd));
	CHECK(is_isomorphic_unrooted(snd, fth));
	CHECK(!is_isomorphic_unrooted(fst.tree, fth));
}

TEST_CASE("is_isomorphic_rooted (complex)", "[validation]") {
	auto t = parse_new_nwk("((((s2,s4),((s13,s1),s7)),s3),s5);");
	CHECK(is_isomorphic_rooted(t.tree,
	                           parse_nwk("((((s2,s4),((s13,s1),s7)),s3),s5);", t.indices)));
	CHECK(is_isomorphic_rooted(t.tree,
	                           parse_nwk("((((s4,s2),((s1,s13),s7)),s3),s5);", t.indices)));
	CHECK(is_isomorphic_rooted(t.tree,
	                           parse_nwk("(s5,(s3,((s2,s4),(s7,(s13,s1)))));", t.indices)));
	CHECK(!is_isomorphic_rooted(t.tree,
	                            parse_nwk("((((s2,s1),((s13,s4),s7)),s3),s5);", t.indices)));
}

} // namespace tests
} // namespace terraces
