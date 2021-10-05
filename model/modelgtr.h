//
//modelgtr.h
//Created by James Barbetti on 06-Jul-2021
//
#include "modelmarkov.h"
#include "utils/eigendecomposition.h"

class ModelGTR : public ModelMarkov {
    public:
        typedef ModelMarkov super;
        ModelGTR(PhyloTree *tree, bool count_rates);
        ~ModelGTR() = default;

        virtual void        saveCheckpoint() override;
        virtual void        restoreCheckpoint() override;
        virtual std::string getName()       const override;
        virtual std::string getNameParams() const override;
        void                getNameParamsFreq(std::ostream &retname) const;
        void                init(const char *model_name, const std::string& model_params,
                                 StateFreqType freq, const std::string& freq_params,
                                 PhyloTree* report_to_tree) override;
        virtual void        writeInfo(ostream &out) override;
        virtual void        getRateMatrix(double *rate_mat) const override;
        virtual void        setRateMatrix(double* rate_mat) override;
        virtual void        setStateFrequency(double* freq) override;
        virtual void        getQMatrix(double *q_mat) const override;
        virtual int         getNDim() const override;
        virtual int         getNDimFreq() const override;
        virtual bool        scaleStateFreq() override;
        virtual void        setVariables(double *variables) override;
        virtual bool        getVariables(const double *variables) override;
        virtual double      targetFunk(double x[]) override;
        virtual bool        isUnstableParameters() const override;
        virtual void        setBounds(double *lower_bound, double *upper_bound, 
	                                  bool *bound_check) override;
        virtual double      optimizeParameters(double gradient_epsilon,
                                               PhyloTree* report_to_tree) override;
        virtual void        decomposeRateMatrix() override;
        virtual void        manuallyComputeEigenvectors();
        virtual void        readRates(istream &in) override;
        virtual void        readRates(std::string str) override;
        virtual void        readStateFreq(istream &in, PhyloTree* report_to_tree) override;
        virtual void        readStateFreq(string str, PhyloTree* report_to_tree) override;
        virtual void        readParameters(const char* file_name, 
                                           bool        adapt_tree_ignored,
							               PhyloTree*  report_to_tree);
        virtual void        freeMem() override;

        virtual void    computeTransMatrix(double time, 
                                   double *trans_matrix, int mixture=0) 
                                   const override;
        void    computeTransMatrixFreq(double time, 
                                       double* trans_matrix) const;
        double  computeTrans(double time, int state1, int state2) const override;
        double  computeTrans(double time, int state1, int state2, 
                             double &derv1, double &derv2) const override;
        virtual void computeTransDerv(double  time, double* trans_matrix, 
	                             double* trans_derv1, double* trans_derv2,
                                 int mixture = 0) const override;
        void    computeTransDervFreq(double time, double rate_val, 
                                     double* trans_matrix, double* trans_derv1, 
				 					 double* trans_derv2) const;

    protected:
        bool       half_matrix;
        int        highest_freq_state;
};