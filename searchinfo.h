/*
 * searchinfo.h
 *
 *  Created on: Nov 18, 2014
 *      Author: Tung
 */

#ifndef SEARCHINFO_H_
#define SEARCHINFO_H_

class SearchInfo {
public:
	SearchInfo();
	void init(bool nni5, double initPS, int maxNNISteps, bool reduction, bool speedNNI);
	virtual ~SearchInfo();

	/* Getters and Setters function */
	double getCurPs() const;
	void setCurPs(double curPs);
	double getInitPs() const;
	void setInitPs(double initPs);
	int getMaxNniSteps() const;
	void setMaxNniSteps(int maxNniSteps);
	bool isNni5() const;
	void setNni5(bool nni5);
	bool isReduction() const;
	void setReduction(bool reduction);
	bool isSpeedNni() const;
	void setSpeedNni(bool speedNni);

	bool nniSkipped;
private:
	bool speedNNI;
	bool reduction;
	bool nni5;
	int maxNNISteps;
	double initPS; // initial perturbation strength
	double curPS; // current perturbation strength

};

#endif /* SEARCHINFO_H_ */
