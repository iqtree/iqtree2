//
// possibleplacement.cpp
// Implementation of the PossiblePlacement class.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "possibleplacement.h"


PossiblePlacement::PossiblePlacement()
                    : parsimony_score(0), score(0)
                    , lenToNewTaxon(-1)
                    , lenToNode1(0), lenToNode2(0) {
}
PossiblePlacement& PossiblePlacement::operator= (const PossiblePlacement & rhs) = default;
bool PossiblePlacement::operator < (const PossiblePlacement& rhs) const {
    return score < rhs.score;
}
bool PossiblePlacement::operator <= (const PossiblePlacement& rhs) const {
    return score <= rhs.score;
}
void PossiblePlacement::setTargetBranch(TargetBranchRange* targetRange, size_t index) {
    target_branch = TargetBranchRef(targetRange, index);
}
void PossiblePlacement::setTargetBranch(TargetBranchRef& branch_ref) {
    target_branch = branch_ref;
}
bool PossiblePlacement::canStillUse() const {
    return !target_branch.isUsedUp();
}
TargetBranch* PossiblePlacement::getTarget() const {
    return target_branch.getTarget();
}
size_t PossiblePlacement::getTargetIndex() const {
    return target_branch.getTargetIndex();
}
void PossiblePlacement::forget() {
    target_branch = TargetBranchRef(nullptr, 0);
}
